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
* File: UQSolver1D.cpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description:  This header file is for the 1-dimensional UQ Solver.
*  This class implements all the numerical methods for solving PDEs, as well as
*  all Monte Carlo methods in one space dimension.
*
* (see also UQSolver1D.hpp)
* ==============================================================================
*/
#include "UQSolver1D.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <boost/timer/progress_display.hpp>
#include <chrono>
#include <random>

#include <iomanip> // @TODO: Get rid of this later

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConstants.hpp"
#include "WriteToFile.hpp"
#include "MemoryAllocations.hpp"
#include "CMode1D.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;
using namespace memory;

/*==============================================================================
    Constructor
==============================================================================*/
/** Main constructor 
 *  @param aModes: vector containing object for each mode 
 *  @param aLambdaMatrix: 2D matrix of transition rates (lambda) 
 *  @param aLambdaEstimates: 3D matrix of bounds on transition rates (lambda) 
 *  @param aFilename: prefix to be used when writing data to file 
 *  @param aFishing:  boolean used to set whether solving fishing example or not */
UQSolver1D::UQSolver1D(vector<CMode1D> aModes,
                       shared_ptr<array2D_t<double>> aLambdaMatrix,
                       shared_ptr<array3D_t<double>> aLambdaEstimates,
                       const string aFilename, const bool aFishing) {
  /** Read in arguments */
  fModes = aModes;
  fM = fModes.size();
  fLambdaMatrix = aLambdaMatrix;
  fLambdaEstimates = aLambdaEstimates;
  fFilename = aFilename;
  fFishing = aFishing;

  /** Set grid parameters */
  fNx   = fModes[0].getGridSizeX();
  fNs   = fModes[0].getGridSizeS();
  fDx   = fModes[0].getDx();
  fDs   = fModes[0].getDs();
  fMinX = fModes[0].getMinX();
  fMaxX = fModes[0].getMaxX();
  fMaxS = fModes[0].getMaxS();

  /** Define stencil to test control values */
  const vector<double> control_values = {-1.0, 0.0, 1.0};
  const int n_controls = control_values.size();

  /** Check CFL-like condition */
  double C_min = LARGE_NUMBER;
  double f_max = -LARGE_NUMBER;
  for (int i = 0; i < fM; i++) {
    for (int k = 0; k < fNx; k++) {
      const double x = fModes[i].xGridToPhysical(k);
      for (int l = 0; l < n_controls; l++) {
        C_min = min(C_min, fModes[i].getCost(x,control_values[l]));
        f_max = max(f_max, abs(fModes[i].getVelocity(x,control_values[l])));
      }
    }
  }
  assert(C_min > 0);
  assert(fDs / C_min <= fDx / f_max);

  /** Write grid and step sizes to File */
  write_config_to_file();

  return;
}

/*==============================================================================
    Monte Carlo main method
==============================================================================*/
/** Compute aNTrials Monte Carlo trials starting from aStartingLocation.
  *  @param aInitialPosition: vector containing starting location
  *  @param aInitialMode:     starting mode
  *  @param aNtrials:         number of Monte Carlo Trials
  *  @param aControlStrategy: choice of control policy to use in MC trials.
                  choices are "CDF", "EV", and "none" */
void UQSolver1D::monteCarlo(const std::vector<double> aInitialPosition,
                            const int aInitialMode, const int aNTrials,
                            const string aControlStrategy) {
  if (aNTrials > 0) {
    /** Set up progress bar and timer for MC trials */
    auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Running Monte Carlo Trials: ";
    boost::timer::progress_display loading_bar(aNTrials);

    /** Get control strategy "deadline" */
    double deadline = fMaxS;
    if (aControlStrategy == "EV") {
      deadline = 0.0;
    } else if (aControlStrategy.rfind("deadline", 0) == 0) {
      /** Get substring containing deadline and convert to double */
      const int index = aControlStrategy.rfind("_");
      if (index >= 0) {
        const string deadline_str = aControlStrategy.substr(index+1,aControlStrategy.length());
        deadline = stod(deadline_str);
      } else {
        std::cerr << "Invalid deadline for Monte Carlo trials" << '\n';
      }
    }

    /** Main computational loop */
    vector<double> mc_results = vector<double>(aNTrials);
    for (int n_trial = 0; n_trial < aNTrials; n_trial++) {
      mc_results[n_trial] = monteCarloTrial(aInitialPosition[0],aInitialMode,deadline);
      ++loading_bar;
    }

    /** Stop timer */
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;

    /** Set up progress bar and timer for MC CDF calculations*/
    auto t3 = std::chrono::high_resolution_clock::now();
    std::cout << "Calculating CDF and EV from Monte Carlo Trials: ";
    boost::timer::progress_display loading_bar2(fNs);

    /** Calculate CDF */
    std::vector<double> cdf_vector = std::vector<double>(fNs);
    for (int j = 0; j < fNs; j++) {
      const double s = fModes[0].sGridToPhysical(j);
      double sum1 = 0.0;

      for (int trial = 0; trial < aNTrials; trial++) {
        if (mc_results[trial] <= s) {
          sum1 += 1.0;
        }
      }

      const double cdf = sum1 / aNTrials;
      cdf_vector[j] = cdf;
      ++loading_bar2;
    }

    /** Calculate EV */
    double sum2 = 0.0;
    for (int trial = 0; trial < aNTrials; trial++) {
      sum2 += mc_results[trial];
    }
    const double avg = sum2 / aNTrials;

    /** Write to file */
    std::vector<int> config = {aNTrials, fMaxSteps};
    io::writeVectorToFile<int>(fFilename + "_MonteCarloConfig", config);
    io::writeVectorToFile<double>(fFilename + "_MC_pt_cdf_" + aControlStrategy, cdf_vector);
    const std::vector<double> pt_avg = {avg};
    io::writeVectorToFile<double>(fFilename + "_MC_pt_avg_" + aControlStrategy, pt_avg);

    /** Stop timer */
    auto t4 = std::chrono::high_resolution_clock::now();
    cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t4-t3).count() << " milliseconds." << endl;
  }
}

/*==============================================================================
    Monte Carlo helper functions
==============================================================================*/
/** Individual Monte Carlo trial.
 *  Given starting position and mode, returns accumulated cost. 
 *   @param aInitialPosition: starting location
 *   @param aInitialMode:     starting mode
 *   @param aDeadline: control strategy will maximize probability of arriving
 *                     before aDeadline, then switch to minimizing EV */
double UQSolver1D::monteCarloTrial(const double aInitialPosition,
                                   const int aInitialMode,
                                   const double aDeadline) const {
  assert(fLambdaMatrix != NULL);

  /** Step size */
  const double dt = fDs;

  /** Create random number generator */
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  std::random_device rd;
  std::mt19937 random_engine{rd()};

  /** Define stencil to test control values */
  const vector<double> control_values = {-1.0, 1.0};
  const int n_controls = control_values.size();

  /** Starting position, time, and mode */
  double x = aInitialPosition;
  int i    = aInitialMode;

  /** Elapsed time and accumulated cost */
  int steps = 0;
  double J = 0;

  /** Main loop */
  while (x > fMinX && (x < fMaxX || (fFishing && x == fMaxX)) && steps < fMaxSteps) {
    double opt_cdf = -LARGE_NUMBER;
    double opt_ev = LARGE_NUMBER;
    double opt_control = 0.0;

    /** Compute remaining budget until deadline */
    const double s = max(aDeadline - J,0.0);

    /** Test control values */
    for (int l = 0; l < n_controls; l++) {
      const double a = control_values[l];
      const vector<double> new_values = mc_update(x,s,i,control_values[l],dt);
      const double new_cdf = new_values[0];
      const double new_ev  = new_values[1];
      if (new_cdf > opt_cdf) {
        opt_cdf = new_cdf;
        opt_ev  = new_ev;
        opt_control = a;
      } else if (new_cdf == opt_cdf) {
        /** Tie-breakers based on expected value */
        if (new_ev < opt_ev) {
          opt_ev = new_ev;
          opt_control = a;
        }
      }
    }

    /** Get velocity and cost */
    const double f = fModes[i].getVelocity(x,opt_control);
    const double C = fModes[i].getCost(x,opt_control);

    /** Increment step */
    steps++;
    J += C*dt;
    x += f*dt;

    /** Calculate transition probabilities */
    double random_number = uniform(random_engine);
    for (int j = 0; j < fM; j++) {
      const double p = transitionProbability(i,j,dt);
      if (random_number <= p) {
        i = j;
        break;
      } else {
        random_number -= p;
      }
    }
  }

  return J;
}

/** Update a single CDF value at specified physical coordinates
  *  @param aX: physical x-coordinate
  *  @param aS: physical s-coordinate
  *  @param aMode: current mode
  *  @param aControl: control value
  *  @param aTau: timestep
  *  output: (cdf_value, ev_value) */
vector<double> UQSolver1D::mc_update(const double aX, const double aS, 
                                     const int aMode, const double aControl,
                                     const double aTau) const {
  /** Get velocity and costs */
  const double f = fModes[aMode].getVelocity(aX,aControl);
  const double C = fModes[aMode].getCost(aX,aControl);

  double new_cdf = 0.0;
  double new_ev  = C * aTau;

  /** Calculate movement */
  const double x1 = aX + f * aTau;
  const double s1 = max(aS - C * aTau, 0.0);

  if (x1 <= fMinX || (!fFishing && x1 >= fMaxX)) {
    /** Made it to the boundary */
    new_cdf = 1.0;
  } else {
    for (int j = 0; j < fM; j++) {
      new_cdf += transitionProbability(aMode,j,aTau) * fModes[j].interpolateCDF(x1,s1);
      new_ev  += transitionProbability(aMode,j,aTau) * fModes[j].interpolateEV(x1,s1);
    }
  }
  return {new_cdf, new_ev};
}
/*==============================================================================
    CDF Calculations
==============================================================================*/
/** Compute the CDF of the piecewise deterministic Markov process in 1D.
 *  @param aUseMinCost: flag for deciding whether to initialize with
 *                      minimum cost computations
 *  @param aControlStrategy: choice of control policy to use .
                  choices are "CDF", "EV", and "none"
                    CDF: choose controls to optimize CDF
                    EV: use EV-optimal controls
                    none: use uniformly zero controls */
void UQSolver1D::computeCDF(const bool aUseMinCost, const string aControlStrategy) {
  assert(fLambdaMatrix != NULL);

  /** If necessary, initialize with minimum cost parameters */
  const bool bounds = false;
  if (aUseMinCost) {
    compute_min_cost(bounds,aControlStrategy);
  }

  /** Set up progress bar and timer */
  auto t1 = std::chrono::high_resolution_clock::now();
  cout << "Computing optimal CDF using PDE numerics: ";
  boost::timer::progress_display loading_bar((fNx-2) * (fNs-1) * fM);

  /** Allocate arrays */
  allocate_variables({"CDF"});

  /** Define stencil to test control values */
  const vector<double> control_values = {-1.0, 1.0};
  const int n_controls = control_values.size();

  /** Initial and boundary conditions */
  for (int i = 0; i < fM; i++) {
    for (int n = 1; n < fNs; n++) {
      for (int k = 0; k < fNx; k++) {
        if (k == 0) {
          fModes[i].setCDF(k,n,1.0);
          fModes[i].setCDF_Control(k,n,-1.0);
        } else {
          fModes[i].setCDF(k, n, 0.0);
          fModes[i].setCDF_Control(k,n,0.0);
        }
      }
    }
    for (int k = 0; k < fNx; k++) {
      fModes[i].setCDF_Control(k,0,fModes[i].getEV_Control(k));
    }
  }

  /** For non-fishing examples, set right boundary */
  if (!fFishing) {
    for (int i = 0; i < fM; i++) {
      for (int n = 0; n < fNs; n++) {
        fModes[i].setCDF(fNx-1,n,1.0);
        fModes[i].setCDF_Control(fNx-1,n,1.0);
      }
    }
  }

  /** Maximum index of particle position */
  int kMax;
  if (fFishing){
    kMax = fNx;
  } else {
    kMax = fNx - 1;
  }

  /** Main Computational Loop */
  /** k = position, i = starting mode, j = ending mode, n = "time", l = control value */
  for (int n = 0; n < fNs - 1; n++) {
    for (int k = 1; k < kMax; k++) {
      for (int i = 0; i < fM; i++) {
        /** Calculate index for minimum cost scenario */
        int n1;
        if (aUseMinCost) {
          n1 = ceil(fMinCost[k] / fDs);
        } else {
          n1 = 0;
        }

        /** Update CDF */
        if (n1 > n + 1 && aUseMinCost) {
          /** Hopeless region H */
          fModes[i].setCDF(k,n+1,0.0);
          fModes[i].setExpectedValue(k,n+1,fModes[i].getExpectedValue(k,0));
          fModes[i].setCDF_Control(k,n+1,fModes[i].getEV_Control(k));
        } else if (n1 == n + 1 && aUseMinCost) {
          /** Boundary of hopeless region H */
          fModes[i].setCDF(k,n+1,(*fMinCostProb)[k][i]);
          fModes[i].setExpectedValue(k,n+1,fModes[i].getExpectedValue(k,0));
          if (aControlStrategy == "CDF") {
            fModes[i].setCDF_Control(k,n+1,fMinCost_Controls[k]);
          } else {
            fModes[i].setCDF_Control(k,n+1,fModes[i].getEV_Control(k));
          }
        } else {
          if (aControlStrategy == "CDF") {
            /** Initialize variables */
            double opt_cdf = -LARGE_NUMBER;
            double opt_ev  = LARGE_NUMBER;
            double opt_control = 0.0;

            /** Test control values */
            for (int l = 0; l < n_controls; l++) {
              const vector<double> new_values = update_cdf(k,n+1,i,control_values[l]);
              const double new_cdf = new_values[0];
              const double new_ev  = new_values[1];

              if (new_cdf > opt_cdf) {
                opt_cdf = new_cdf;
                opt_ev  = new_ev;
                opt_control = control_values[l];
              } else if (new_cdf == opt_cdf) {
                /** Tie-breakers based on expected value */
                if (new_ev < opt_ev) {
                  /** Check expected value */
                  opt_ev = new_ev;
                  opt_control = control_values[l];
                }
              }
            }
            fModes[i].setCDF(k,n+1,opt_cdf);
            fModes[i].setExpectedValue(k,n+1,opt_ev);
            fModes[i].setCDF_Control(k,n+1,opt_control);
          } else {
            /** Use fixed control */
            double a = 0.0;
            if (aControlStrategy == "EV") {
              a = fModes[i].getEV_Control(k);
            }
            const vector<double> new_values = update_cdf(k,n+1,i,a);
            const double new_cdf = new_values[0];
            const double new_ev  = new_values[1];
            fModes[i].setCDF(k,n+1,new_cdf);
            fModes[i].setExpectedValue(k,n+1,new_ev);
            fModes[i].setCDF_Control(k,n+1,a);
          }
        }
        ++loading_bar;
      }
    }
  }
  /** Stop timer */
  auto t2 = std::chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;

  /** Write to file */
  write_modes_to_file({"CDF", "EV"});
}

/** Update a single CDF value at specified logical coordinates
  *  @param aK: logical x-coordinate
  *  @param aN: logical s-coordinate
  *  @param aMode: current mode
  *  @param aControl: control value
  *  output: (cdf_value, ev_value) */
vector<double> UQSolver1D::update_cdf(const int aK, const int aN, 
                                      const int aMode, 
                                      const double aControl) const {
  /** Get local coordinates */
  const double x = fModes[aMode].xGridToPhysical(aK);

  /** Get velocity and costs */
  const double f = fModes[aMode].getVelocity(x,aControl);
  const double C = fModes[aMode].getCost(x,aControl);

  /** Calculate timestep */
  const double tau = fDs / C;

  double new_cdf = 0.0;
  double new_ev = C * tau;

  /** Calculate movement */
  const double x1 = x + f * tau;
  const double s1 = fModes[aMode].sGridToPhysical(aN-1);

  if (x1 <= fMinX || (!fFishing && (x1 >= fMaxX))) {
    /** Made it to the boundary */
    new_cdf = 1.0;
  } else {
    for (int j = 0; j < fM; j++) {
      new_cdf += transitionProbability(aMode,j,tau) * fModes[j].interpolateCDF(x1,s1);
      new_ev  += transitionProbability(aMode,j,tau) * fModes[j].interpolateEV(x1,s1);
    }
  }
  return {new_cdf, new_ev};
}

/*==============================================================================
    Bounds on CDF
==============================================================================*/
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
void UQSolver1D::computeCDFBounds(const bool aUseMinCost,
                                  const string aControlStrategy,
                                  const bool aNumBound) {
  assert(!fFishing);
  assert(fLambdaEstimates != NULL);

  const bool bounds = true;
  if (aUseMinCost) {
    compute_min_cost(bounds,aControlStrategy);
  }

  /** Set up progress bar and timer */
  auto t1 = std::chrono::high_resolution_clock::now();
  cout << "Computing bounds on CDF using PDE numerics: ";
  boost::timer::progress_display loading_bar((fNx-2) * (fNs-1) * fM);

  /** Allocate arrays */
  allocate_variables({"CDFBounds"});

  /** Define stencil to test control values */
  const vector<double> control_values = {-1.0, 1.0};
  const int n_controls = control_values.size();

  /** Initial and boundary conditions */
  for (int i = 0; i < fM; i++) {
    for (int n = 0; n < fNs; n++) {
      for(int k = 0; k < fNx; k++) {
        if (k == 0) {
          fModes[i].setCDFUpperBound(k,n,1.0);
          fModes[i].setCDFLowerBound(k,n,1.0);
          fModes[i].setCDF_Control(k,n,-1.0);
        } else {
          fModes[i].setCDFUpperBound(k,n,0.0);
          fModes[i].setCDFLowerBound(k,n,0.0);
        }
      }
    }
  }

  /** For non-fishing examples, set right boundary */
  if (!fFishing) {
    for (int i = 0; i < fM; i++) {
      for (int n = 0; n < fNs; n++) {
        fModes[i].setCDFUpperBound(fNx-1,n,1.0);
        fModes[i].setCDFLowerBound(fNx-1,n,1.0);
        fModes[i].setCDF_Control(fNx-1,n,1.0);
      }
    }
  }

  /** Maximum index of particle position */
  int kMax;
  if (fFishing){
    kMax = fNx;
  } else {
    kMax = fNx - 1;
  }

  /** Main Computational Loop */
  /** k = position, i = starting mode, j = ending mode, n = "time", l = control value */
  for (int n = 0; n < fNs - 1; n++) {
    for (int k = 1; k < kMax; k++) {
      for (int i = 0; i < fM; i++) {
        /** Calculate index for minimum cost scenario */
        int n1;
        if (aUseMinCost) {
          n1 = ceil(fMinCost[k] / fDs);
        } else {
          n1 = 0;
        }

        /** Update CDF */
        if (n1 > n + 1 && aUseMinCost) {
          /** Hopeless region H */
          fModes[i].setCDFUpperBound(k,n+1,0.0);
          fModes[i].setCDFLowerBound(k,n+1,0.0);
        } else if (n1 == n + 1 && aUseMinCost) {
          /** Boundary of hopeless region H */
          fModes[i].setCDFUpperBound(k,n+1,(*fMinCostProbUpper)[k][i]);
          fModes[i].setCDFLowerBound(k,n+1,(*fMinCostProbLower)[k][i]);
          fModes[i].setCDF_Control(k,n+1,fMinCost_Controls[k]);
        } else {
          if (aControlStrategy == "CDF") {
            /** Initialize variables */
            double opt_cdf_upper = -LARGE_NUMBER;
            double opt_cdf_lower = -LARGE_NUMBER;
            double opt_control = 0.0;

            /** Test control values */
            /** Note: current implementation only stores controls for minimizing
             *        lower bound, i.e. controlling worst case scenario */
            for (int l = 0; l < n_controls; l++) {
              /** To calculate optimized (minimized) lower bound on CDF*/
              const double new_cdf_lower = update_cdf_bounds(k,n+1,i,control_values[l],0, aNumBound);
              if (new_cdf_lower > opt_cdf_lower) {
                opt_cdf_lower = new_cdf_lower;
                opt_control = control_values[l];
              }
              /** To calculate optimized (minimized) upper bound on CDF */
              const double new_cdf_upper = update_cdf_bounds(k,n+1,i,control_values[l],1, aNumBound);
              if (new_cdf_upper > opt_cdf_upper) {
                opt_cdf_upper = new_cdf_upper;
              }
            }
            fModes[i].setCDFUpperBound(k,n+1,opt_cdf_upper);
            fModes[i].setCDFLowerBound(k,n+1,opt_cdf_lower);
            fModes[i].setCDF_Control(k,n+1,opt_control);
          } else {
            /** Use fixed control */
            double a = 0.0;
            if (aControlStrategy == "EV") {
              a = fModes[i].getCDF_Control(k,n+1);
            }
            const double new_cdf_lower = update_cdf_bounds(k,n+1,i,a,0, aNumBound);
            const double new_cdf_upper = update_cdf_bounds(k,n+1,i,a,1, aNumBound);
            fModes[i].setCDFUpperBound(k,n+1,new_cdf_upper);
            fModes[i].setCDFLowerBound(k,n+1,new_cdf_lower);
          }
        }
        ++loading_bar;
      }
    }
  }
  /** Stop timer */
  auto t2 = std::chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;

  /** Write to file */
  write_modes_to_file({"CDFBounds"});
}

/** Update bounds for a single CDF value
  *  @param aK: logical x-coordinate
  *  @param aN: logical s-coordinate
  *  @param aMode: current mode
  *  @param aControl: control value
  *  @param aB: determines whether updating lower (0) or upper (1) bound
  *  @param aNumBound: flag for deciding whether to compute numerical or conservative bounds
  *                    True: use numerical bounds
  *                    False: use conservative analytical bounds */
double UQSolver1D::update_cdf_bounds(const int aK, const int aN, const int aMode,
                                     const double aControl, const int aB, 
                                     const bool aNumBound) const {
  double new_cdf = 0.0;
  double cdf_if_stay = 0.0;
  double cdf_if_switch;
  const double x = fModes[aMode].xGridToPhysical(aK);
  const double f = fModes[aMode].getVelocity(x,aControl);
  const double C = fModes[aMode].getCost(x,aControl);

  const double tau = fDs / C;

  /** Calculate movement */
  const double x1 = x + f*tau;
  const double s1 = fModes[aMode].sGridToPhysical(aN-1);
  if (x1 <= fMinX || (!fFishing && x1 >= fMaxX)) {
    /** Made it to the boundary */
    new_cdf = 1.0;
  } else {
    /** Under a linearity assumption on the relationship between lambda and
      * the transition probabilities, we need not check every combination
      * of upper and lower bounds, but can just compare the bounds on the
      * cdfs in the various modes to the bounds on the cdf assuming we stay. */
    if(aB == 0) {
      cdf_if_stay = fModes[aMode].interpolateCDFLowerBound(x1,s1);
    } else if(aB == 1) {
      cdf_if_stay = fModes[aMode].interpolateCDFUpperBound(x1,s1);
    }

    /** Sets cdf_if_stay to be estimate if no mode change occurs */
    new_cdf = cdf_if_stay;

    for (int j = 0; j < fM; j++) {
      if (j != aMode) {
        if (aB == 0) {
          /** Lower bound */
          if (aNumBound) {
            cdf_if_switch = fModes[j].interpolateCDFLowerBound(x1,s1);
          } else {
            cdf_if_switch = 0;
          }
          if (cdf_if_switch > cdf_if_stay) {
            new_cdf += transitionProbabilityEstimate(aMode,j,tau,0)*(cdf_if_switch-cdf_if_stay);
          } else {
            new_cdf += transitionProbabilityEstimate(aMode,j,tau,1)*(cdf_if_switch-cdf_if_stay);
          }
        } else if (aB == 1) {
          /** Upper bound */
          if (aNumBound) {
            cdf_if_switch = fModes[j].interpolateCDFUpperBound(x1,s1);
          } else {
            cdf_if_switch = 1;
          }
          if (cdf_if_switch > cdf_if_stay) {
            new_cdf += transitionProbabilityEstimate(aMode,j,tau,1)*(cdf_if_switch-cdf_if_stay);
          } else {
            new_cdf += transitionProbabilityEstimate(aMode,j,tau,0)*(cdf_if_switch-cdf_if_stay);
          }
        }
      }
    }
  }
  return new_cdf;
}

/*==============================================================================
    Expected Value computations
==============================================================================*/
/** Compute expectation-optimal policy and corresponding expected value
 *  using locking sweeping method.
 *  @param aUseFixedControl: flag for deciding whether to use fixed controls
 *                             True:  use uniform zero controls
 *                             False: chooses controls to optimize EV 
 *                                    and places these controls in fModes */
void UQSolver1D::computeExpectedValue(const bool aUseFixedControl) {
  assert(!fFishing);
  assert(fLambdaMatrix != NULL);

  /** Set up timer */
  auto t1 = std::chrono::high_resolution_clock::now();
  cout << "Computing expected value via sweeping method ... " << flush;

  /** Initialize array for holding average case results */
  for (int k = 0; k < fNx; k++) {
    for (int i = 0; i < fM; i++) {
      if (k == 0) {
        fModes[i].setExpectedValue(k,0,0.0);
        fModes[i].setEV_Control(k,-1.0);
      } else {
        fModes[i].setExpectedValue(k,0,LARGE_NUMBER);
      }
    }
  }

  /** For non-fishing examples, set right boundary */
  if (!fFishing) {
    for (int i = 0; i < fM; i++) {
      fModes[i].setExpectedValue(fNx-1,0,0.0);
      fModes[i].setEV_Control(fNx-1,1.0);
    }
  }

  /** Initialize locks */
  array2D_t<bool> active = allocateArray2D<bool>(fNx,fM);
  for (int k = 0; k < fNx; k++) {
    for (int i = 0; i < fM; i++) {
      if (k <= 1 || k >= fNx-2) {
        active[k][i] = true;
      } else {
        active[k][i] = false;
      }
    }
  }

  /** Convergence condition */
  constexpr double eps = 1.0e-6;

  /** Sweep directions */
  vector<int> sweep_inc  = {1,-1};
  vector<int> sweep_init = {1,fNx-2};

  /** Main Computational Loop */
  bool converged = false;
  while (!converged) {
    double max_change = 0.0;
    /** Enumerate sweep directions */
    for (int sweep = 0; sweep < 2; sweep++) {
      /** Iterate over gridpoints */
      for (int k = sweep_init[sweep]; k > 0 && k < fNx - 1; k += sweep_inc[sweep]) {
        /** Iterate over modes */
        for (int i = 0; i < fM; i++) {
          /** Check if gridpoint is active */
          if (active[k][i]) {
            /** Set gridpoint as inactive */
            active[k][i] = false;

            /** Try to update value */
            const double old_value = fModes[i].getExpectedValue(k,0);
            const bool was_updated = update_ev(k,i,aUseFixedControl);

            if (was_updated) {
              /** Calculate change */
              const double new_value = fModes[i].getExpectedValue(k,0);
              const double change = old_value - new_value;
              max_change = max(change, max_change);

              /** Update locks for neighbors */
              if (k > 1) {
                for (int j = 0; j < fM; ++j) {
                  active[k-1][j] = true;
                }
              }
              if (k < fNx - 2) {
                for (int j = 0; j < fM; ++j) {
                  active[k+1][j] = true;
                }
              }
            }
          }
        }
      }
    }

    /** Check convergence */
    if (max_change < eps) {
      converged = true;
    }
  }

  /** Stop timer */
  cout << "Done." << endl;
  auto t2 = std::chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;
}

/** Helper method for Expected value calculations.
  * Update expected value at gridpoint (k,i) */
bool UQSolver1D::update_ev(const int aK, const int aI, 
                           const bool aUseFixedControl) {
  bool was_updated = false;

  /** Define stencil to test control values */
  vector<double> control_values;
  if (aUseFixedControl) {
    control_values = {0.0};
  } else {
    control_values = {-1.0, 0.0, 1.0};
  }
  const int n_controls = control_values.size();

  /** Get local coordinates */
  const double x = fModes[aI].xGridToPhysical(aK);

  double update = LARGE_NUMBER;
  double opt_control = 0.0;
  /** Iterate over control values */
  for (int l = 0; l < n_controls; l++) {
    const double a = control_values[l];

    /** Get cost and speed functions */
    const double C = fModes[aI].getCost(x,a);
    const double f = fModes[aI].getVelocity(x,a);

    /** Check zero speed case */
    double temp;
    if (f == 0) {
      double numerator = C;
      double denominator = 0.0;
      for (int i = 0; i < fM; i++) {
        if (i != aI) {
          numerator += (*fLambdaMatrix)[aI][i] * fModes[i].getExpectedValue(aK,0);
          denominator += (*fLambdaMatrix)[aI][i];
        }
      }
      temp = numerator / denominator;
    } else {
      /** Calculate semi-Lagrangian update */
      const double tau = fDx / abs(f);
      int k1;
      if (f > 0) {
        k1 = aK + 1;
      } else {
        k1 = aK - 1;
      }
      temp = C * tau;
      for (int i = 0; i < fM; i++) {
        temp += transitionProbability(aI,i,tau) * fModes[i].getExpectedValue(k1,0);
      }
    }
    if (temp < update) {
      opt_control = a;
      update = temp;
    }
  }

  /** Only update if new value actually decreases EV */
  const double old_value = fModes[aI].getExpectedValue(aK,0);
  if (update < old_value) {
    fModes[aI].setExpectedValue(aK, 0, update);
    fModes[aI].setEV_Control(aK, opt_control);
    was_updated = true;
  }
  return was_updated;
}

/*==============================================================================
    Minimum Cost computations
==============================================================================*/
/** Compute minimum cost and minimum cost probabilities using Djikstra-like method
  *  @param aBounds - notes whether or not we are calculating bounds on CDF 
  *  @param aControlStrategy: choice of control policy to use.
                  choices are "CDF", "EV", and "none" */
void UQSolver1D::compute_min_cost(const bool aBounds, 
                                  const string aControlStrategy) {
  /** Set up timer */
  auto t1 = chrono::high_resolution_clock::now();
  cout << "Computing minimum costs and probabilities ... " << flush;

  /** Allocate vector/arrays for minimum cost computations */
  fMinCost          = vector<double>(fNx);
  fMinCost_Controls = vector<double>(fNx);
  fMinCostProbLower = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fM));
  fMinCostProbUpper = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fM));
  fMinCostProb      = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fM));

  /** Initialize heap and arrays for Djikstra's */
  heap_t heap;
  vector<status_t> status(fNx);
  vector<handle_t> heap_pointers(fNx);
  for (int k = 0; k < fNx; k++) {
    if (!fFishing && (k == fNx - 1)) {
      fMinCost[k] = 0.0;
      heap.push(HeapGP(k,0.0));
      status[k] = ACCEPTED;
      for (int i = 0; i < fM; i++) {
        (*fMinCostProbUpper)[k][i] = 1.0;
        (*fMinCostProbLower)[k][i] = 1.0;
        (*fMinCostProb)[k][i] = 1.0;
      }
    } else if (k == 0) {
      fMinCost[k] = 0.0;
      heap.push(HeapGP(k,0.0));
      status[k] = ACCEPTED;
      for (int i = 0; i < fM; i++) {
        (*fMinCostProbUpper)[k][i] = 1.0;
        (*fMinCostProbLower)[k][i] = 1.0;
        (*fMinCostProb)[k][i] = 1.0;
      }
    } else {
      fMinCost[k] = LARGE_NUMBER;
      status[k] = FAR;
      for (int i = 0; i < fM; i++) {
        (*fMinCostProbUpper)[k][i] = 0.0;
        (*fMinCostProbLower)[k][i] = 0.0;
        (*fMinCostProb)[k][i] = 0.0;
      }
    }
  }

  const vector<int> stencil = {1,-1};
  /** Perform Djikstra's algorithm to compute s^0 */
  while (!heap.empty()) {
    /** Grab smallest element of heap */
    HeapGP current_gp = heap.top();
    heap.pop();
    const int k = current_gp.k;
    status[k] = ACCEPTED;

    /** Iterate over neighbors */
    for (int index = 0; index < 2; index++) {
      /** Get coordinates of neighbor */
      const int k1 = k + stencil[index];

      const bool in_bounds = (k1 >= 0 && k1 <= fNx-1);

      /** If accepted nothing to do, otherwise, try to update */
      if (in_bounds && (status[k1] != ACCEPTED)) {
        /** Update gridpoint */
        bool was_updated = update_min_cost(k1, aBounds, aControlStrategy);
        /** If it was already considered and was updated then update heap */
        if (status[k1] == CONSIDERED) {
          if (was_updated) {
            /** We want to use increase here because the boost implementation is a max heap */
            heap.increase(heap_pointers[k1], HeapGP(k1, fMinCost[k1]));
          }
        } else {
          /** Else add to heap */
          status[k1] = CONSIDERED;
          heap_pointers[k1] = heap.push(HeapGP(k1, fMinCost[k1]));
        }
      }
    }
  }
  /** Stop timer */
  auto t2 = chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;

  /** Write results to file */
  io::writeVectorToFile<double>(fFilename + "_MinCost", fMinCost);
  io::writeToFile2D<double>(fFilename + "_MinCostProb", *fMinCostProb);
  io::writeVectorToFile<double>(fFilename + "_MinCostControls", fMinCost_Controls);
}

/** Helper method for Minimum cost calculations.
  * Update expected value at gridpoint (k) */
bool UQSolver1D::update_min_cost(const int aK, const bool aBounds, 
                                 const string aControlStrategy) {
  /** Set flag for whether or not point was updated */
  bool was_updated = false;

  /** Calculate coordinates */
  const double x = fModes[0].xGridToPhysical(aK);

  /** Define stencil to test control values */
  std::vector<double> control_values;
  if (aControlStrategy == "CDF") {
    control_values = {-1.0, 1.0};
  } else {
    control_values = {0.0};
  }
  const int n_controls = control_values.size();

  double opt_value = fMinCost[aK];

  /** Iterate over all modes */
  for (int i = 0; i < fM; i++) {
    /** Minimize costs over all controls */
    for (int l = 0; l < n_controls; l++) {
      double a = control_values[l];
      if (aControlStrategy == "EV") {
        a = fModes[i].getEV_Control(aK);
      }
      const double C = fModes[i].getCost(x,a);
      const double f = fModes[i].getVelocity(x,a);

      double new_value;
      if (f == 0) {
        new_value = LARGE_NUMBER;
      } else {
        /** Find indices of neighbor */
        int k1;
        if (f >= 0) {
          k1 = aK + 1;
        } else {
          k1 = aK - 1;
        }

        /** Update value */
        const double dt = fDx / abs(f);
        new_value = fMinCost[k1] + C*dt;

        /** Compare new update to existing value */
        if (new_value <= opt_value) {
          /** Update minimum cost probabilites for bounds */
          if (aBounds) {
            (*fMinCostProbUpper)[aK][i] = (*fMinCostProbUpper)[k1][i];
            (*fMinCostProbLower)[aK][i] = (*fMinCostProbLower)[k1][i];
            for (int j = 0; j < fM; j++) {
              const double pL = transitionProbabilityEstimate(i,j,dt, 0);
              const double pU = transitionProbabilityEstimate(i,j,dt, 1);
              const double lower_diff = (*fMinCostProbLower)[k1][j] - (*fMinCostProbLower)[k1][i];
              const double upper_diff = (*fMinCostProbUpper)[k1][j] - (*fMinCostProbUpper)[k1][i];
              (*fMinCostProbLower)[aK][i] += min(pL * lower_diff, pU * lower_diff);
              (*fMinCostProbUpper)[aK][i] += max(pL * upper_diff, pU * upper_diff);
              if (j != i && new_value < opt_value) {
                (*fMinCostProbLower)[aK][j] = 0.0;
                (*fMinCostProbUpper)[aK][j] = 0.0;
              }
            }
          } else {
            /** Update minimum cost and minimum cost probabilities */
            (*fMinCostProb)[aK][i] = (*fMinCostProb)[k1][i];
            for (int j = 0; j < fM; j++) {
              const double p = transitionProbability(i,j,dt);
              (*fMinCostProb)[aK][i] += p * ((*fMinCostProb)[k1][j] - (*fMinCostProb)[k1][i]);
              if (j != i && new_value < opt_value) {
                (*fMinCostProb)[aK][j] = 0.0;
              }
            }
          }

          if (new_value < opt_value) {
            opt_value = new_value;
          }
          fMinCost_Controls[aK] = a;

          was_updated = true;
        }
      }
    }
  }

  /** Update minimum cost */
  fMinCost[aK] = opt_value;

  return was_updated;
}

/*==============================================================================
    Write to file
==============================================================================*/
/** Write grid sizes, step sizes, and lambda values to file */
void UQSolver1D::write_config_to_file() const {
  /** Write lambda matrix to file */
  if (fLambdaMatrix != NULL) {
    io::writeToFile2D<double>(fFilename + "_lambda", *fLambdaMatrix);
  }

  if (fLambdaEstimates != NULL) {
    io::writeToFile3D<double>(fFilename + "_lambdaEstimates", *fLambdaEstimates);
  }

  /** Write grid sizes to file */
  std::vector<int> gridsize = {fNx, fNs, fM};
  io::writeVectorToFile<int>(fFilename + "_Gridsizes", gridsize);

  /** Write step sizes to file */
  std::vector<double> config = {fDx, fDs, fMinX, fMaxX, fMaxS};
  io::writeVectorToFile<double>(fFilename + "_Stepsizes", config);
}

/** Write CDF and Expected Values to file */
void UQSolver1D::write_modes_to_file(const vector<string> aVariables) const {
  /** Write grids to file */
  for (int j = 0; j < fM; ++j) {
    fModes[j].writeModeToFile(aVariables, fFilename + std::to_string(j));
  }
}

/*==============================================================================
    Allocating arrays
==============================================================================*/
/** An I/O member which writes the grid functions for each mode to file.
 * @param aVariables vector of variables to write to file.
 *      Choices are "CDF" and "CDFBounds" */
void UQSolver1D::allocate_variables(const std::vector<std::string> aVariables) {
  for (int j = 0; j < fM; ++j) {
    fModes[j].allocateVariables(aVariables);
  }
}