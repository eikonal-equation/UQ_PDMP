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
* File: UQSolver2D.cpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This file contains the numerical methods for computing CDFs
*   in two space dimensions.
*
* (see also UQSolver2D.hpp)
* ==============================================================================
*/
#include "UQSolver2D.hpp"

/** ------ Libraries ---------------------------------------------------------*/
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>
#include <boost/timer/progress_display.hpp>
#include <chrono>
#include <random>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConstants.hpp"
#include "WriteToFile.hpp"
#include "MemoryAllocations.hpp"
#include "CMode2D.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;
using namespace std::placeholders;
using namespace memory;

/*==============================================================================
    Constructor
==============================================================================*/
/** Main constructor 
 *  @param aModes: vector containing object for each mode 
 *  @param aLambdaMatrix: 2D matrix of transition rates (lambda) 
 *  @param aFilename: prefix to be used when writing data to file */
UQSolver2D::UQSolver2D(vector<CMode2D> aModes,
                       shared_ptr<array2D_t<double>> aLambdaMatrix,
                       const string aFilename) {
  /** Read in arguments */
  fModes = aModes;
  fM = fModes.size();
  fLambdaMatrix = aLambdaMatrix;
  fFilename = aFilename;

  /** Read in grid parameters */
  fNx   = fModes[0].getGridSizeX();
  fNy   = fModes[0].getGridSizeY();
  fNs   = fModes[0].getGridSizeS();
  fDx   = fModes[0].getDx();
  fDy   = fModes[0].getDy();
  fDs   = fModes[0].getDs();
  fMinX = fModes[0].getMinX();
  fMinY = fModes[0].getMinY();
  fMaxX = fModes[0].getMaxX();
  fMaxY = fModes[0].getMaxY();
  fMaxS = fModes[0].getMaxS();

  /** Set number of controls for manually testing CFL-like condition */
  constexpr int n_controls = 120;
  vector<double> theta_values = {};
  for (int i = 0; i < n_controls; i++) {
    theta_values.push_back(2 * PI * i / n_controls);
  }

  /** Check CFL-like condition */
  double C_min = LARGE_NUMBER;
  double f_max = -LARGE_NUMBER;
  for (int k = 0; k < fM; k++) {
    for (int i = 0; i < fNx; i++) {
      for (int j = 0; j < fNy; j++) {
        for (int l = 0; l < n_controls; l++) {
          const double a = theta_values[l];
          const double x = fModes[k].xGridToPhysical(i);
          const double y = fModes[k].yGridToPhysical(j);
          const double K = fModes[k].getCost(x, y, a);
          const double f_x = fModes[k].getVelocityX(x, y, a);
          const double f_y = fModes[k].getVelocityY(x, y, a);
          const double f = max(abs(f_x),abs(f_y));
          C_min = min(C_min, K);
          f_max = max(f_max, f);
        }
      }
    }
  }
  assert(C_min > 0);
  assert(fDs / C_min <= min(fDx,fDy) / f_max);

  /** Write grid and step sizes to File */
  write_config_to_file();

  return;
}

/*==============================================================================
    Monte Carlo main method
==============================================================================*/
/** Compute aNTrials Monte Carlo trials starting from aInitialPosition.
  *  @param aInitialPosition: vector containing starting location
  *  @param aInitialMode:     starting mode
  *  @param aNtrials:         number of Monte Carlo Trials
  *  @param aControlStrategy: choice of control policy to use in MC trials.
                  choices are "CDF", "EV", and "none" */
void UQSolver2D::monteCarlo(const std::vector<double> aInitialPosition,
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
      mc_results[n_trial] = monteCarloTrial(aInitialPosition,aInitialMode,deadline);
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
 *   @param aInitialPosition: vector containing starting location
 *   @param aInitialMode:     starting mode
 *   @param aDeadline:        control strategy will maximize probability of arriving
 *                            before aDeadline, then switch to minimizing EV */
double UQSolver2D::monteCarloTrial(const std::vector<double> aInitialPosition,
                                   const int aInitialMode, 
                                   const double aDeadline) const {
  /** Set step size */
  const double dt = fDs;

  /** Create random number generator */
  std::uniform_real_distribution<double> uniform(0.0,1.0);
  std::random_device rd;
  std::mt19937 random_engine{rd()};

  /** Starting position and mode */
  assert(aInitialPosition.size() >= 2);
  double x = aInitialPosition[0];
  double y = aInitialPosition[1];
  int i    = aInitialMode;

  /** Elapsed time and accumulated cost */
  int steps = 0;
  double J = 0;

  /** Main loop */
  while (x > fMinX && x < fMaxX && y > fMinY && y < fMaxY && steps < fMaxSteps) {
    /** Compute remaining budget until deadline */
    const double s = max(aDeadline - J,0.0);

    /** Optimize over control values */
    const vector<double> new_values = mc_optimize(x,y,s,i,dt);
    const double a_opt = new_values[2];

    /** Get velocity and cost */
    const double f_x = fModes[i].getVelocityX(x, y, a_opt);
    const double f_y = fModes[i].getVelocityY(x, y, a_opt);
    const double C = fModes[i].getCost(x, y, a_opt);

    /** Increment step */
    steps++;
    J += C*dt;
    x += f_x*dt;
    y += f_y*dt;

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

/** Choose control to optimize CDF value
  *  @param aX: physical x-coordinate
  *  @param aY: physical y-coordinate
  *  @param aS: physical s-coordinate
  *  @param aMode:    current mode
  *  @param aTau:     timestep
  *  output: optimal (cdf,ev,control) */
vector<double> UQSolver2D::mc_optimize(const double aX, const double aY,
                                       const double aS, const int aMode,
                                       const double aTau) const {
  double opt_cdf = -LARGE_NUMBER;
  double opt_ev  = LARGE_NUMBER;
  double opt_grid_control;

  /** Define stencil for testing control values */
  constexpr int n_controls = 8;
  vector<double> theta_values = {};
  for (int i = 0; i < n_controls; i++) {
    theta_values.push_back(2 * PI * i / n_controls);
  }

  /** Test grid of control values */
  for (int l = 0; l < n_controls; l++) {
    const double a = theta_values[l];
    const vector<double> new_values = mc_update(aX,aY,aS,aMode,a,aTau);
    const double new_cdf = new_values[0];
    const double new_ev  = new_values[1];
    if (new_cdf > opt_cdf) {
      opt_cdf = new_cdf;
      opt_ev  = new_ev;
      opt_grid_control = a;
    } else if (new_cdf == opt_cdf) {
      /** Tiebreakers */
      if (new_ev < opt_ev) {
        /** First check expected value */
        opt_ev = new_ev;
        opt_grid_control = a;
      } 
    }
  }

  /** Optimize using golden section search */
  const double lower_bound = opt_grid_control - 2*PI/n_controls;
  const double upper_bound = opt_grid_control + 2*PI/n_controls;
  const std::function<std::vector<double>(double)> objective_function = std::bind(&UQSolver2D::mc_update, this, aX, aY, aS, aMode, _1, aTau);
  const double a_opt = cdf_golden_section_search(objective_function, lower_bound, upper_bound);
  const vector<double> opt_values = mc_update(aX,aY,aS,aMode,a_opt,aTau);
  return {opt_values[0], opt_values[1], a_opt};
}

/** Update a single CDF value at physical coordinates
  *  @param aX: physical x-coordinate
  *  @param aY: physical y-coordinate
  *  @param aS: physical s-coordinate
  *  @param aMode:    current mode
  *  @param aControl: control value
  *  @param aTau:     timestep
  *  output: (cdf_value, ev_value) */
vector<double> UQSolver2D::mc_update(const double aX, const double aY,
                                     const double aS, const int aMode,
                                     const double aControl,
                                     const double aTau) const {
  /** Get velocity and costs */
  const double f_x = fModes[aMode].getVelocityX(aX,aY,aControl);
  const double f_y = fModes[aMode].getVelocityY(aX,aY,aControl);
  const double C   = fModes[aMode].getCost(aX,aY,aControl);

  double new_cdf = 0.0;
  double new_ev  = C * aTau;

  /** Calculate movement */
  const double x1 = aX + f_x * aTau;
  const double y1 = aY + f_y * aTau;
  const double s1 = max(aS - C * aTau, 0.0);

  if (x1 <= fMinX || x1 >= fMaxX || y1 <= fMinY || y1 >= fMaxY) {
    /** Made it to the edge */
    new_cdf = 1.0;
  } else {
    for (int j = 0; j < fM; j++) {
      if (aMode != j) {
        /** Transition */
        new_cdf += transitionProbability(aMode,j,aTau) * fModes[j].interpolateCDF(x1,y1,s1);
        new_ev  += transitionProbability(aMode,j,aTau) * fModes[j].interpolateEV(x1,y1,s1);
      } else {
        /** Stay in same mode */
        new_cdf += transitionProbability(j,j,aTau) * fModes[j].interpolateCDF(x1,y1,s1);
        new_ev  += transitionProbability(j,j,aTau) * fModes[j].interpolateEV(x1,y1,s1);
      }
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
void UQSolver2D::computeCDF(const bool aUseMinCost, 
                            const string aControlStrategy) {
  /** If necessary, initialize with best case parameters */
  if (aUseMinCost) {
    compute_min_cost(aControlStrategy);
  }

  /** Set up progress bar and timer */
  auto t1 = chrono::high_resolution_clock::now();
  cout << "Computing CDF using PDE numerics: " << flush;
  boost::timer::progress_display loading_bar((fNx-2) * (fNy-2) * (fNs-1) * fM);

  /** Initial and boundary conditions */
  for (int k = 0; k < fM; k++) {
    for (int n = 0; n < fNs; n++) {
      for (int i = 0; i < fNx; i++) {
        for (int j = 0; j < fNy; j++) {
          if (i == 0) {
            /** Left boundary */
            fModes[k].setCDF(i,j,n,1.0);
            fModes[k].setCDF_Control(i,j,n,PI);
          } else if (i == fNx - 1) {
            /** Right boundary */
            fModes[k].setCDF(i,j,n,1.0);
            fModes[k].setCDF_Control(i,j,n,0.0);
          } else if (j == 0) {
            /** Lower boundary */
            fModes[k].setCDF(i,j,n,1.0);
            fModes[k].setCDF_Control(i,j,n,3*PI/2);
          } else if (j == fNy - 1) {
            /** Upper boundary */
            fModes[k].setCDF(i,j,n,1.0);
            fModes[k].setCDF_Control(i,j,n,PI/2);
          } else {
            fModes[k].setCDF(i, j, n, 0.0);
          }
        }
      }
    }
    for (int i = 0; i < fNx; i++) {
      for (int j = 0; j < fNy; j++) {
        fModes[k].setCDF_Control(i,j,0,fModes[k].getEV_Control(i,j));
      }
    }
  }

  /** Main Computational Loop */
  /** k = starting mode, i = x position, j = y position, k0 = ending mode, n = "time" */
  for (int n = 0; n < fNs - 1; n++) {
    for (int i = 1; i < fNx - 1; i++) {
      for (int j = 1; j < fNy - 1; j++) {
        for (int k = 0; k < fM; k++) {
          /** Calculate index for minimum cost scenario */
          int n1;
          if (aUseMinCost) {
            n1 = ceil((*fMinCost)[i][j] / fDs);
          } else {
            n1 = 0;
          }

          /** Update CDF */
          if (n1 > n + 1 && aUseMinCost) {
            /** Hopeless region H */
            fModes[k].setCDF(i,j,n+1,0.0);
            fModes[k].setExpectedValue(i,j,n+1,fModes[k].getExpectedValue(i,j,0));
            const double a = fModes[k].getEV_Control(i,j);
            fModes[k].setCDF_Control(i,j,n+1,a);
          } else if (n1 == n + 1 && aUseMinCost) {
            /** Boundary of hopeless region H */
            fModes[k].setCDF(i,j,n+1,(*fMinCostProb)[i][j][k]);
            fModes[k].setExpectedValue(i,j,n+1,fModes[k].getExpectedValue(i,j,0));
            if (aControlStrategy == "CDF") {
              const double a = (*fMinCost_Controls)[i][j];
              fModes[k].setCDF_Control(i,j,n+1,a);
            } else {
              const double a = fModes[k].getEV_Control(i,j);
              fModes[k].setCDF_Control(i,j,n+1,a);
            }
          } else {
            if (aControlStrategy == "CDF") {
              const vector<double> opt_values = optimize_cdf(i,j,n+1,k);

              /** Set new values */
              fModes[k].setCDF(i,j,n+1,opt_values[0]);
              fModes[k].setExpectedValue(i,j,n+1,opt_values[1]);
              fModes[k].setCDF_Control(i,j,n+1,opt_values[2]);
            } else {
              /** Use fixed control */
              double a = 0.0;
              if (aControlStrategy == "EV") {
                a = fModes[k].getEV_Control(i,j);
              }

              const vector<double> new_values = update_cdf(i,j,n+1,k,a);
              const double new_cdf = new_values[0];
              const double new_ev  = new_values[1];
              fModes[k].setCDF(i,j,n+1,new_cdf);
              fModes[k].setExpectedValue(i,j,n+1,new_ev);
              fModes[k].setCDF_Control(i,j,n+1,a);
            }
          }
          ++loading_bar;
        }
      }
    }
  }
  /** Stop timer */
  auto t2 = chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;

  /** Write to file */
  write_modes_to_file({"CDF", "EV"});
}

/** Choose control to optimize CDF value
  *  @param aI: logical x-coordinate
  *  @param aJ: logical y-coordinate
  *  @param aN: logical s-coordinate
  *  @param aMode:    current mode
  *  output: optimal (cdf,ev,control) */
vector<double> UQSolver2D::optimize_cdf(const int aI, const int aJ,
                                        const int aN, const int aMode) const {
  double opt_cdf = -LARGE_NUMBER;
  double opt_ev  = LARGE_NUMBER;
  double opt_grid_control;

  /** Define stencil for testing control values */
  constexpr int n_controls = 8;
  vector<double> theta_values = {};
  for (int i = 0; i < n_controls; i++) {
    theta_values.push_back(2 * PI * i / n_controls);
  }

  /** Test coarse grid of control values */
  for (int l = 0; l < n_controls; l++) {
    const double a = theta_values[l];
    const vector<double> new_values = update_cdf(aI,aJ,aN,aMode,a);
    const double new_cdf = new_values[0];
    const double new_ev  = new_values[1];
    if (new_cdf > opt_cdf) {
      opt_cdf = new_cdf;
      opt_ev  = new_ev;
      opt_grid_control = a;
    } else if (new_cdf == opt_cdf) {
      /** Tiebreakers */
      if (new_ev < opt_ev) {
        /** First check expected value */
        opt_ev = new_ev;
        opt_grid_control = a;
      } 
    }
  }

  /** Optimize using golden section search */
  const double lower_bound = opt_grid_control - 2*PI/n_controls;
  const double upper_bound = opt_grid_control + 2*PI/n_controls;
  const std::function<std::vector<double>(double)> objective_function = std::bind(&UQSolver2D::update_cdf, this, aI, aJ, aN, aMode, _1);
  const double a_opt = cdf_golden_section_search(objective_function, lower_bound, upper_bound);
  const vector<double> opt_values = update_cdf(aI,aJ,aN,aMode,a_opt);
  return {opt_values[0], opt_values[1], a_opt};
}

/** Recursive implementation of golden section search (for CDF maximization)
 *    This is designed to include tiebreakers,
 *  @param f:    function to be maximized 
 *              should take in a double and output vector containing two doubles
 *              this method will then maximize first output
 *              while using secound output as a tiebreaker
 *  @param a,b:  endpoints of interval on which to use golden section search
 *  returns: maximizer x */
double UQSolver2D::cdf_golden_section_search(const std::function<std::vector<double>(double)> f,
                                             double a, double b) const {
  double c = b - (b - a) / GOLDEN_RATIO;
  double d = a + (b - a) / GOLDEN_RATIO;
  while (abs(c - d) > fTolGSS) {
    const vector<double> f_c = f(c);
    const vector<double> f_d = f(d);
    if (f_c[0] > f_d[0] || (f_c[0] == f_d[0] && f_c[1] < f_d[1])) {
      b = d;
    } else {
      a = c;
    }
    c = b - (b - a) / GOLDEN_RATIO;
    d = a + (b - a) / GOLDEN_RATIO;
  }
  return (b + a) / 2;
}

/** Update a single CDF value
  *  @param aI: logical x-coordinate
  *  @param aJ: logical y-coordinate
  *  @param aN: logical s-coordinate
  *  @param aMode:    current mode
  *  @param aControl: control value
  *  output: (cdf_value, ev_value) */
vector<double> UQSolver2D::update_cdf(const int aI, const int aJ,
                                      const int aN, const int aMode,
                                      const double aControl) const {
  /** Get local coordinates */
  const double x = fModes[aMode].xGridToPhysical(aI);
  const double y = fModes[aMode].yGridToPhysical(aJ);

  /** Get velocity and costs */
  const double f_x = fModes[aMode].getVelocityX(x,y,aControl);
  const double f_y = fModes[aMode].getVelocityY(x,y,aControl);
  const double C   = fModes[aMode].getCost(x,y,aControl);

  const double tau = fDs / C;

  double new_cdf = 0.0;
  double new_ev  = C * tau;

  /** Calculate movement */
  const double x1 = x + f_x*tau;
  const double y1 = y + f_y*tau;
  const double s1 = fModes[aMode].sGridToPhysical(aN-1);

  if (x1 <= fMinX || x1 >= fMaxX || y1 <= fMinY || y1 >= fMaxY) {
    /** Made it to the edge */
    new_cdf = 1.0;
  } else {
    for (int j = 0; j < fM; j++) {
      new_cdf += transitionProbability(aMode,j,tau) * fModes[j].interpolateCDF(x1,y1,s1);
      new_ev  += transitionProbability(aMode,j,tau) * fModes[j].interpolateEV(x1,y1,s1);
    }
  }
  return {new_cdf, new_ev};
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
void UQSolver2D::computeExpectedValue(const bool aUseFixedControl) {
  assert(fDx == fDy);

  /** Set up timer */
  auto t1 = std::chrono::high_resolution_clock::now();
  cout << "Computing expected value via sweeping method ... " << flush;

  /** Initialize array for holding average case results */
  for (int i = 0; i < fNx; i++) {
    for (int j = 0; j < fNy; j++) {
      for (int k = 0; k < fM; k++) {
        if (i == 0) {
          fModes[k].setExpectedValue(i,j,0,0.0);
          fModes[k].setEV_Control(i,j,PI);
        } else if (i == fNx-1) {
          fModes[k].setExpectedValue(i,j,0,0.0);
          fModes[k].setEV_Control(i,j,0.0);
        } else if (j == 0) {
          fModes[k].setExpectedValue(i,j,0,0.0);
          fModes[k].setEV_Control(i,j,3*PI/2);
        } else if (j == fNy-1) {
          fModes[k].setExpectedValue(i,j,0,0.0);
          fModes[k].setEV_Control(i,j,PI/2);
        } else {
          fModes[k].setExpectedValue(i,j,0,LARGE_NUMBER);
        }
      }
    }
  }

  /** Initialize locks */
  array3D_t<bool> active = allocateArray3D<bool>(fNx,fNy,fM);
  for (int i = 0; i < fNx; i++) {
    for (int j = 0; j < fNy; j++) {
      for (int k = 0; k < fM; k++) {
        if (i <= 1 || i >= fNx -2 || j <= 1 || j >= fNy-2) {
          active[i][j][k] = true;
        } else {
          active[i][j][k] = false;
        }
      }
    }
  }

  /** Convergence condition */
  constexpr double eps = 1.0e-6;

  /** Sweep directions */
  vector<int> sweep_inc_x  = {1,-1,-1,1};
  vector<int> sweep_inc_y  = {1,1,-1,-1};
  vector<int> sweep_init_x = {1,fNx-2,fNx-2,1};
  vector<int> sweep_init_y = {1,1,fNy-2,fNy-2};

  /** Main Computational Loop */
  bool converged = false;
  while (!converged) {
    double max_change = 0.0;
    /** Enumerate sweep directions */
    for (int sweep = 0; sweep < 4; sweep++) {
      /** Iterate over gridpoints */
      for (int i = sweep_init_x[sweep]; i > 0 && i < fNx - 1; i += sweep_inc_x[sweep]) {
        for (int j = sweep_init_y[sweep]; j > 0 && j < fNy - 1; j += sweep_inc_y[sweep]) {
          /** Iterate over modes */
          for (int k = 0; k < fM; ++k) {
            /** Check if gridpoint is active */
            if (active[i][j][k]) {
              /** Set gridpoint as inactive */
              active[i][j][k] = false;

              /** Try to update value */
              const double old_value = fModes[k].getExpectedValue(i,j,0);
              vector<double> new_values;
              if (aUseFixedControl) {
                const double a = 0.0;
                new_values = {update_ev(i,j,k,a), a};
              } else {
                new_values = optimize_ev(i,j,k);
              }

              if (new_values[0] < old_value) {
                /** Set new value function and controls */
                fModes[k].setExpectedValue(i,j,0,new_values[0]);
                fModes[k].setEV_Control(i,j,new_values[1]);

                /** Calculate change */
                const double change = old_value - new_values[0];
                max_change = max(change, max_change);

                /** Update locks for neighbors */
                for (int k1 = 0; k1 < fM; ++k1) {
                  if (i > 1) {
                    active[i-1][j][k1] = true;
                  }
                  if (i < fNx - 2) {
                    active[i+1][j][k1] = true;
                  }
                  if (j > 1) {
                    active[i][j-1][k1] = true;
                  }
                  if (j < fNy - 2) {
                    active[i][j+1][k1] = true;
                  }
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

/** Recursive implementation of golden section search (for EV maximization)
 *  @param f:    function to be maximized
 *               f should be a function with one input and one output
 *  @param a,b:  endpoints of interval
 *  returns: maximizer x */
double UQSolver2D::ev_golden_section_search(const std::function<double(double)> f,
                                            double a, double b) const {
  double c = b - (b - a) / GOLDEN_RATIO;
  double d = a + (b - a) / GOLDEN_RATIO;
  while (abs(c - d) > fTolGSS) {
    if (f(c) < f(d)) {
      b = d;
    } else {
      a = c;
    }
    c = b - (b - a) / GOLDEN_RATIO;
    d = a + (b - a) / GOLDEN_RATIO;
  }
  return (b + a) / 2;
}

vector<double> UQSolver2D::optimize_ev(const int aI, const int aJ, const int aMode) const {
  double opt_grid_ev = LARGE_NUMBER;
  double opt_grid_control;

  /** Define coarse stencil for testing control values */
  constexpr int n_controls = 8;
  vector<double> theta_values = {};
  for (int i = 0; i < n_controls; i++) {
    theta_values.push_back(2 * PI * i / n_controls);
  }

  /** Test coarse grid of control values */
  for (int l = 0; l < n_controls; l++) {
    const double a = theta_values[l];
    const double new_ev = update_ev(aI,aJ,aMode,a);
    if (new_ev < opt_grid_ev) {
      opt_grid_ev = new_ev;
      opt_grid_control = a;
    }
  }

  /** Optimize using golden section search */
  const double lower_bound = opt_grid_control - 2*PI/n_controls;
  const double upper_bound = opt_grid_control + 2*PI/n_controls;
  const std::function<double(double)> objective_function = std::bind(&UQSolver2D::update_ev, this, aI, aJ, aMode, _1);
  const double a_opt = ev_golden_section_search(objective_function, lower_bound, upper_bound);
  const double opt_ev = update_ev(aI,aJ,aMode,a_opt);
  return {opt_ev, a_opt};
}

/** Compute EV update at gridpoint (aI,aJ) and mode aMode using aControl */
double UQSolver2D::update_ev(const int aI, const int aJ, const int aMode, 
                             const double aControl) const {
  /** Get local coordinates */
  const double x = fModes[aMode].xGridToPhysical(aI);
  const double y = fModes[aMode].yGridToPhysical(aJ);

  /** Get cost and speed functions */
  const double C   = fModes[aMode].getCost(x,y,aControl);
  const double f_x = fModes[aMode].getVelocityX(x,y,aControl);
  const double f_y = fModes[aMode].getVelocityY(x,y,aControl);

  /** Check zero speed case */
  if (f_x == 0 && f_y == 0) {
    return LARGE_NUMBER;
  }

  /** Find indices of vertical and horizontal neighbors */
  int i1, j1;
  if (f_x >= 0) {
    i1 = aI + 1;
  } else {
    i1 = aI - 1;
  }
  if (f_y >= 0) {
    j1 = aJ + 1;
  } else {
    j1 = aJ - 1;
  }

  /** Calculate semi-Lagrangian update */
  const double eta = abs(f_y)/(abs(f_x) + abs(f_y));
  const double dt = fDx / (abs(f_x) + abs(f_y));

  double update = C * dt;
  for (int k = 0; k < fM; ++k) {
    const double u_x = fModes[k].getExpectedValue(i1,aJ,0);
    const double u_y = fModes[k].getExpectedValue(aI,j1,0);
    const double u_interp = eta * u_y + (1 - eta) * u_x;
    update += transitionProbability(aMode,k,dt) * u_interp;
  }

  return update;
}

/*==============================================================================
    Minimum cost computations
==============================================================================*/
/** Compute minimum cost and associated probabilities using Djikstra-like method
  *  @param aControlStrategy: choice of control policy to use.
                  choices are "CDF", "EV", and "none" */
void UQSolver2D::compute_min_cost(const string aControlStrategy) {
  assert(fDx == fDy);

  /** Set up timer */
  auto t1 = chrono::high_resolution_clock::now();
  cout << "Computing minimum costs and probabilities ... " << flush;

  /** Allocate vector/arrays for minimum cost computations */
  fMinCost          = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));
  fMinCostProb      = make_shared<array3D_t<double>>(allocateArray3D<double>(fNx,fNy,fM));
  fMinCost_Controls = make_shared<array2D_t<double>>(allocateArray2D<double>(fNx,fNy));

  /** Initialize heap and arrays for Djikstra's */
  heap_t heap;
  array2D_t<status_t> status        = allocateArray2D<status_t>(fNx,fNy);
  array2D_t<handle_t> heap_pointers = allocateArray2D<handle_t>(fNx,fNy);
  for (int i = 0; i < fNx; i++) {
    for (int j = 0; j < fNy; j++) {
      if (i == 0 || i == fNx - 1 || j == 0 || j == fNy - 1) {
        (*fMinCost)[i][j] = 0.0;
        heap.push(HeapGP(i,j,0.0));
        status[i][j] = ACCEPTED;
        for (int k = 0; k < fM; k++) {
          (*fMinCostProb)[i][j][k] = 1.0;
        }
      } else {
        (*fMinCost)[i][j] = LARGE_NUMBER;
        status[i][j] = FAR;
        for (int k = 0; k < fM; k++) {
          (*fMinCostProb)[i][j][k] = 0.0;
        }
      }
    }
  }

  /** Vectors for storing 5pt stencil */
  const vector<double> stencil_x = {1,-1,0,0};
  const vector<double> stencil_y = {0,0,-1,-1};

  /** Perform Djikstra's algorithm to compute s^0 */
  while (!heap.empty()) {
    /** Grab smallest element of heap */
    HeapGP current_gp = heap.top();
    heap.pop();
    const int i = current_gp.i;
    const int j = current_gp.j;
    status[i][j] = ACCEPTED;

    /** Update neighbors */
    for (int index = 0; index < 4; index++) {
      /** Get coordinates of neighbor */
      const int i1 = i + stencil_x[index];
      const int j1 = j + stencil_y[index];

      const bool in_bounds = (i1 >= 0 && i1 <= fNx-1 && j1 >= 0 && j1 <= fNy-1);

      /** If accepted nothing to do, otherwise, try to update */
      if (in_bounds && (status[i1][j1] != ACCEPTED)) {
        /** Update gridpoint */
        bool was_updated = update_min_cost(i1,j1,aControlStrategy);

        /** If it was already considered and was updated then update heap */
        if (status[i1][j1] == CONSIDERED) {
          if (was_updated) {
            /** We want to use increase here because the boost implementation is a max heap */
            heap.increase(heap_pointers[i1][j1], HeapGP(i1, j1, (*fMinCost)[i1][j1]));
          }
        } else {
          /** Else add to heap */
          status[i1][j1] = CONSIDERED;
          heap_pointers[i1][j1] = heap.push(HeapGP(i1, j1, (*fMinCost)[i1][j1]));
        }
      }
    }
  }
  /** Stop timer */
  auto t2 = chrono::high_resolution_clock::now();
  cout << "Took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds." << endl;

  /** Write results to file */
  io::writeToFile2D<double>(fFilename + "_MinCost", *fMinCost);
  io::writeToFile3D<double>(fFilename + "_MinCostProb", *fMinCostProb);
  io::writeToFile2D<double>(fFilename + "_MinCost_Controls", *fMinCost_Controls);
}

/** Helper method for Minimum cost calculations.
  * Update EV at gridpoint (aI,aJ) using strategy aControlStrategy
  *    Returns True if value is updated,
  *    Returns False if value is not updated. */
/** Currently assumes dx = dy */
bool UQSolver2D::update_min_cost(const int aI, const int aJ, 
                                 const string aControlStrategy) {
  /** Set flag for whether or not point was updated */
  bool was_updated = false;

  /** Calculate coordinates */
  const double x = fModes[0].xGridToPhysical(aI);
  const double y = fModes[0].yGridToPhysical(aJ);

  /** Define stencil to test control values */
  int n_controls;
  vector<double> theta_values;
  if (aControlStrategy == "CDF") {
    n_controls = 120;
    theta_values = {};
    for (int i = 0; i < n_controls; i++) {
      theta_values.push_back(2 * PI * i / n_controls);
    }
  } else {
    n_controls = 1;
    theta_values = {0.0};
  }

  double opt_value = (*fMinCost)[aI][aJ];

  /** Iterate over all modes */
  for (int k = 0; k < fM; k++) {
    /** Minimize costs over all controls */
    for (int l = 0; l < n_controls; l++) {
      double a = theta_values[l];
      if (aControlStrategy == "EV") {
        a = fModes[k].getEV_Control(aI,aJ);
      }
      const double C   = fModes[k].getCost(x,y,a);
      const double f_x = fModes[k].getVelocityX(x,y,a);
      const double f_y = fModes[k].getVelocityY(x,y,a);

      double new_value;
      if (f_x == 0 && f_y == 0) {
        new_value = LARGE_NUMBER;
      } else {
        /** Find indices of vertical and horizontal neighbors */
        int i1, j1;
        if (f_x >= 0) {
          i1 = aI + 1;
        } else {
          i1 = aI - 1;
        }
        if (f_y >= 0) {
          j1 = aJ + 1;
        } else {
          j1 = aJ - 1;
        }

        /** Interpolate between vertical and horizontal neighbors */
        const double s_x = (*fMinCost)[i1][aJ];
        const double s_y = (*fMinCost)[aI][j1];
        const double eta = abs(f_y)/(abs(f_x) + abs(f_y));
        const double s_interp = eta * s_y + (1 - eta) * s_x;
        const double w_x = (*fMinCostProb)[i1][aJ][k];
        const double w_y = (*fMinCostProb)[aI][j1][k];
        const double w_interp = eta * w_y + (1 - eta) * w_x;

        /** Update value */
        const double dt = fDx / (abs(f_x) + abs(f_y));
        new_value = C*dt + s_interp;

        /** Compare new update to existing value */
        if (new_value <= opt_value) {
          /** Update minimum cost and minimum cost probabilites */
          (*fMinCostProb)[aI][aJ][k] = w_interp;
          for (int k1 = 0; k1 < fM; k1++) {
            const double p = transitionProbability(k,k1,dt);
            const double w_k1_x = (*fMinCostProb)[i1][aJ][k1];
            const double w_k1_y = (*fMinCostProb)[aI][j1][k1];
            const double w_k1_interp = eta * w_k1_y + (1 - eta) * w_k1_x;
            (*fMinCostProb)[aI][aJ][k] += p * (w_k1_interp - w_interp);
            if (k1 != k && new_value < opt_value) {
              (*fMinCostProb)[aI][aJ][k1] = 0.0;
            }
          }

          if (new_value < opt_value) {
            opt_value = new_value;
          }
          (*fMinCost_Controls)[aI][aJ] = a;

          was_updated = true;
        }
      }
    }
  }

  /** Update minimum cost */
  (*fMinCost)[aI][aJ] = opt_value;

  return was_updated;
}

/*==============================================================================
    Write to file
==============================================================================*/
/** Write grid sizes, and all variables to the file at aFilename */
void UQSolver2D::write_config_to_file() const {
  /** Write lambda matrix to file */
  io::writeToFile2D<double>(fFilename + "_lambda", *fLambdaMatrix);

  /** Write grid sizes to file */
  std::vector<int> gridsize = {fNx, fNy, fNs, fM};
  io::writeVectorToFile<int>(fFilename + "_Gridsizes", gridsize);

  /** Write step sizes to file */
  std::vector<double> config = {fDx, fDy, fDs, fMinX, fMaxX, fMinY, fMaxY, fMaxS};
  io::writeVectorToFile<double>(fFilename + "_Stepsizes", config);
}

/** An I/O member which writes the grid functions for each mode to file.
 * @param aVariables vector of variables to write to file.
 *      Choices are "CDF", "EV", "Velocity", and "Cost" */
void UQSolver2D::write_modes_to_file(const vector<string> aVariables) const {
  /** Write CDF grids to file */
  for (int j = 0; j < fM; ++j) {
    fModes[j].writeModeToFile(aVariables, fFilename + std::to_string(j));
  }
}