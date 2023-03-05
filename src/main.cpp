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
* File: main.cpp
*
* Author: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This is the main executable file.
* To run examples 1-7 from the manuscript, run the main method with argument 1-7
*
* ==============================================================================
*/

/** ------ Libraries ---------------------------------------------------------*/
#include <chrono>
#include <iostream>
#include <string>
#include <iomanip>

/** ------ Project-specific header files -------------------------------------*/
#include "SimpleFunctions.hpp"
#include "MemoryAllocations.hpp"
#include "CMode1D.hpp"
#include "CMode2D.hpp"
#include "UQSolver.hpp"
#include "UQSolver1D.hpp"
#include "UQSolver2D.hpp"

/*==============================================================================
        One-dimensional Examples
==============================================================================*/
/** ------ Example 1: --------------------------------------------------------*/
/**   - 1D domain: [0,1]
 *    - Uncontrolled
 *    - Mode 1: Move right with speed 1
 *    - Mode 2: Move left with speed 1  
 *    - Transition rates: lambda_12 = lambda_21 = 2 */
UQSolver1D Basic(const int aNx, const int aNs) {
  /** Define grid parameters */
  constexpr double x_min = 0.0;
  constexpr double x_max = 1.0;
  constexpr double s_max = 2.0;
  constexpr bool fishing = false;

  /** Create grids for each mode */
  constexpr int n_modes = 2;
  CMode1D mode1 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_right);
  CMode1D mode2 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_left);
  std::vector<CMode1D> modes = {mode1, mode2};
  assert(n_modes == modes.size());

  /** Define lambda matrix */
  const double temp_lambda[n_modes][n_modes] = {{0.0,2.0},{2.0,0.0}};
  memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);
  std::shared_ptr<memory::array3D_t<double>> lambda_estimates = NULL;

  std::string filename = "Basic";
  return UQSolver1D(modes, lambda_ptr, lambda_estimates, filename, fishing);
}

/** ------ Example 2: --------------------------------------------------------*/
/**   - 1D domain: [0,1]
 *    - Uncontrolled
 *    - Mode 1: Move right with speed 1
 *    - Mode 2: Move left with speed 0.5
 *    - Transition rates: lambda_12 = lambda_21 = 2 */
UQSolver1D SpeedTest(const int aNx, const int aNs) {
  /** Define grid parameters */
  constexpr double x_min = 0.0;
  constexpr double x_max = 1.0;
  constexpr double s_max = 2.0;
  constexpr bool fishing = false;

  /** Create grids for each mode */
  constexpr int n_modes = 2;
  CMode1D mode1 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::slow_right);
  CMode1D mode2 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_left);
  std::vector<CMode1D> modes = {mode1, mode2};
  assert(n_modes == modes.size());

  /** Define lambda matrix */
  const double temp_lambda[n_modes][n_modes] = {{0.0,2.0},{2.0,0.0}};
  memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);
  std::shared_ptr<memory::array3D_t<double>> lambda_estimates;

  std::string filename = "SpeedTest";
  return UQSolver1D(modes, lambda_ptr, lambda_estimates, filename, fishing);
}

/** ------ Example 4: --------------------------------------------------------*/
/**   - 1D domain: [0,1]
 *    - Uncontrolled
 *    - Mode 1: Move right with speed 1
 *    - Mode 2: Move left with speed 0.5  
 *    Used for computing bounds on CDF, with 1 <= lambda <= 4 */
UQSolver1D CDFBounds(const int aNx, const int aNs) {
  /** Define grid parameters */
  constexpr double x_min = 0.0;
  constexpr double x_max = 1.0;
  constexpr double s_max = 2.0;
  constexpr bool fishing = false;

  /** Create grids for each mode */
  constexpr int n_modes = 2;
  CMode1D mode1 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_right);
  CMode1D mode2 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_left);
  std::vector<CMode1D> modes = {mode1, mode2};
  assert(n_modes == modes.size());

  /** Define lambda estimates matrix */
  const double temp_estimates[2][n_modes][n_modes] = {{{0.0,1.0},{1.0,0.0}},{{0.0,4.0},{4.0,0.0}}};
  memory::array3D_t<double> estimate_matrix = memory::initializeArray3D<double,2,n_modes,n_modes>(temp_estimates);
  std::shared_ptr<memory::array3D_t<double>> lambda_estimates = std::make_shared<memory::array3D_t<double>>(estimate_matrix);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = NULL;

  std::string filename = "CDF_Bounds";

  return UQSolver1D(modes, lambda_ptr, lambda_estimates, filename, fishing);
}

/** ------ Example 4: --------------------------------------------------------*/
/**   - 1D domain: [0,1]
 *    - Uncontrolled
 *    - Mode 1: Move right with speed 1
 *    - Mode 2: Move left with speed 0.5  
 *    Transition rates specified by arguments lambda1 and lambda2 
 *    Used for producing green curves in Figures 7 and 8 */
UQSolver1D AsymmetricLambda(const int aNx, const int aNs,
                            const double lambda1, const double lambda2) {
    /** Define grid parameters */
    constexpr double x_min = 0.0;
    constexpr double x_max = 1.0;
    constexpr double s_max = 2.0;
    constexpr bool fishing = false;

    /** Create grids for each mode */
    constexpr int n_modes = 2;
    CMode1D mode1 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_right);
    CMode1D mode2 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::constant_left);
    std::vector<CMode1D> modes = {mode1, mode2};
    assert(n_modes == modes.size());

    /** Define lambda matrix */
    const double temp_lambda[n_modes][n_modes] = {{0.0,lambda1},{lambda2,0.0}};
    memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
    std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);
    std::shared_ptr<memory::array3D_t<double>> lambda_estimates = NULL;

  	std::ostringstream streamObj1;
  	streamObj1 << std::fixed;
  	streamObj1 << std::setprecision(1);
  	streamObj1 << lambda1;
  	std::string strObj1 = streamObj1.str();

    std::ostringstream streamObj2;
  	streamObj2 << std::fixed;
  	streamObj2 << std::setprecision(1);
  	streamObj2 << lambda2;
  	std::string strObj2 = streamObj2.str();

    std::string filename = "CDF_Asym_Lambda_" + strObj1 + "_" + strObj2 + "_";
    return UQSolver1D(modes, lambda_ptr, lambda_estimates, filename, fishing);
}

/** ------ Example 5: --------------------------------------------------------*/
/**   - 1D domain: [0,1]
 *    - Controlled: -1 <= a <= 1
 *    - Mode 1: Move with velocity  0.5 + a
 *    - Mode 2: Move with velocity -0.5 + a
 *    - Transition rates: lambda_12 = lambda_21 = 2 */
UQSolver1D Controlled(const int aNx, const int aNs,
                      const std::string aControlStrategy) {
  /** Define grid parameters */
  constexpr double x_min = 0.0;
  constexpr double x_max = 1.0;
  constexpr double s_max = 1.0;
  constexpr bool fishing = false;

  /** Create grids for each mode */
  constexpr int n_modes = 2;
  CMode1D mode1 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::control_right);
  CMode1D mode2 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::constant2, dynamics::control_left);
  std::vector<CMode1D> modes = {mode1, mode2};
  assert(n_modes == modes.size());

  /** Define lambda matrix */
  const double temp_lambda[n_modes][n_modes] = {{0.0,2.0},{2.0,0.0}};
  memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);
  std::shared_ptr<memory::array3D_t<double>> lambda_estimates;

  std::string filename;
  if (aControlStrategy == "CDF") {
    filename = "CDF_Controlled";
  } else if (aControlStrategy == "EV") {
    filename = "EV_Controlled";
  } else {
    std::cerr << "Invalid control strategy for Example 5" << std::endl;
    assert(false);
  }

  return UQSolver1D(modes, lambda_ptr, lambda_estimates, filename, fishing);
}

/** ------ Example 7: --------------------------------------------------------*/
/**   - 1D
 *    - @TODO: Fill this out */
UQSolver1D Fishing(const int aNx, const int aNs){
  /** Define grid parameters */
  constexpr double x_min = 1.0;
  constexpr double x_max = 4.0;
  constexpr double s_max = 200.0; // give 3.8 mode time to hit population = 1
  constexpr bool fishing = true;

  /** Create grids for each mode */
  constexpr int n_modes = 3;
  CMode1D mode1 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::harvest_cost, dynamics::fish_low);
  CMode1D mode2 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::harvest_cost, dynamics::fish_mid);
  CMode1D mode3 = CMode1D(aNx, aNs, x_min, x_max, s_max, costs::harvest_cost, dynamics::fish_high);
  std::vector<CMode1D> modes = {mode1, mode2, mode3};
  assert(n_modes == modes.size());

  /** Define lambda matrix */
  const double temp_lambda[n_modes][n_modes] = {{0.0,0.1,0.0},{0.05,0.0,0.05},{0.0,0.1,0.0}};
  memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);
  std::shared_ptr<memory::array3D_t<double>> lambda_estimates = NULL;

  std::string filename = "Fishing";
  return UQSolver1D(modes, lambda_ptr, lambda_estimates, filename, fishing);
}

/*==============================================================================
        Two-Dimensional Examples
==============================================================================*/
/** ------ Example 3: --------------------------------------------------------*/
/**   - 2D domain: [0,1]x[0,1]
 *    - Uncontrolled
 *    - Mode 1: Move with velocity (-1,0)
 *    - Mode 2: Move with velocity (0,1)
 *    - Mode 3: Move with velocity (1,0)
 *    - Mode 4: Move with velocity (0,-1)
 *    - Transition rates: lambda_ij = 1 for all i,j */
UQSolver2D TwoDimensional(const int aNx, const int aNs) {
  /** Define grid parameters */
  constexpr double x_min = 0.0;
  constexpr double x_max = 1.0;
  constexpr double y_min = 0.0;
  constexpr double y_max = 1.0;
  constexpr double s_max = 1.0;

  /** Create grids for each mode */
  constexpr int n_modes = 4;
  CMode2D mode_l = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::negative, dynamics::zero);
  CMode2D mode_u = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::zero, dynamics::positive);
  CMode2D mode_r = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::positive, dynamics::zero);
  CMode2D mode_d = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::zero, dynamics::negative);
  std::vector<CMode2D> modes = {mode_l, mode_u, mode_r, mode_d};
  assert(n_modes == modes.size());

  /** Define lambda matrix */
  const double temp_lambda[n_modes][n_modes] = {{0.0,1.0,1.0,1.0},{1.0,0.0,1.0,1.0},{1.0,1.0,0.0,1.0},{1.0,1.0,1.0,0.0}};
  memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);

  const std::string filename = "TwoDimensional";
  return UQSolver2D(modes, lambda_ptr, filename);
}

/** ------ Example 6: --------------------------------------------------------*/
/**   - 2D domain: [0,1]x[0,1]
 *    - Controlled: 0 <= a <= 2*pi
 *    - Mode 1: Move with velocity (-0.5 + cos(a),  sin(a)      )
 *    - Mode 2: Move with velocity (cos(a)       ,  0.5 + sin(a))
 *    - Mode 3: Move with velocity (0.5  + cos(a),  sin(a)      )
 *    - Mode 4: Move with velocity (cos(a)       , -0.5 + sin(a))
 *    - Transition rates: lambda_ij = 1 for all i,j */
UQSolver2D TwoDimensionalControlled(const int aNx, const int aNs, const std::string aControlStrategy) {
  /** Define grid parameters */
  constexpr double x_min = 0.0;
  constexpr double x_max = 1.0;
  constexpr double y_min = 0.0;
  constexpr double y_max = 1.0;
  constexpr double s_max = 1.0;

  /** Create grids for each mode */
  constexpr int n_modes = 4;
  CMode2D mode_l = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::control_2d_left_x, dynamics::control_2d_left_y);
  CMode2D mode_u = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::control_2d_up_x, dynamics::control_2d_up_y);
  CMode2D mode_r = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::control_2d_right_x, dynamics::control_2d_right_y);
  CMode2D mode_d = CMode2D(aNx, aNx, aNs, x_min, x_max, y_min, y_max, s_max, costs::constant3, dynamics::control_2d_down_x, dynamics::control_2d_down_y);
  std::vector<CMode2D> modes = {mode_l, mode_u, mode_r, mode_d};
  assert(n_modes == modes.size());

  /** Define lambda matrix */
  const double temp_lambda[n_modes][n_modes] = {{0.0,1.0,1.0,1.0},{1.0,0.0,1.0,1.0},{1.0,1.0,0.0,1.0},{1.0,1.0,1.0,0.0}};
  memory::array2D_t<double> lambda_matrix = memory::initializeArray2D<double,n_modes,n_modes>(temp_lambda);
  std::shared_ptr<memory::array2D_t<double>> lambda_ptr = std::make_shared<memory::array2D_t<double>>(lambda_matrix);

  std::string filename;
  if (aControlStrategy == "CDF") {
    filename = "2D_CDF_Controlled";
  } else if (aControlStrategy == "EV") {
    filename = "2D_EV_Controlled";
  } else {
    std::cerr << "Invalid control strategy for Example 6" << std::endl;
    assert(false);
  }

  return UQSolver2D(modes, lambda_ptr, filename);
}

/*==============================================================================
        Main function
==============================================================================*/
int main(int argc, char* argv[]){
  /** Check for argument */
  if (argc < 2) {
    std::cerr << "Error: too few arguments" << std::endl;
    assert(false);
  }

  /** Read in first argument: example number */
  const int example = read_arg(argv[1]);
  if (example < 1 || example > 7) {
    std::cerr << "Invalid example number, must be between 1 and 7" << '\n';
  }

  /** Read in optional grid size arguments */
  int nx, ns, n_trials;
  if (argc > 3) {
    nx = read_arg(argv[2]);
    ns = read_arg(argv[3]);
  } else {
    /** Default grid size values */
    switch(example) {
      case 1:
        nx = 1001;
        ns = 2001;
        n_trials = 100;
        break;
      case 2:
        nx = 1001;
        ns = 2001;
        n_trials = 100;
        break;
      case 3:
        nx = 101;
        ns = 101;
        n_trials = 100;
        break;
      case 4:
        nx = 1001;
        ns = 2001;
        n_trials = 100;
        break;
      case 5:
        nx = 8001;
        ns = 16001;
        n_trials = 100000;
        break;
      case 6:
        nx = 401;
        ns = 601;
        n_trials = 10000;
        break;
      case 7:
        nx = 101;
        ns = 120001;
        n_trials = 10000;
        break;
    }
  }

  /** Initialize variables */
  constexpr bool use_min_cost = true;
  bool fixed_control_ev;
  bool num_bound;
  std::string cdf_control_strategy;
  std::string mc_control_strategy;
  std::shared_ptr<UQSolver> uq_solver;
  std::vector<double> init_location;
  int init_mode;
  std::vector<std::string> deadlines;

  /** Start timer */
  auto t1 = std::chrono::high_resolution_clock::now();
  
  switch(example) {
    case 1:
      /** ------ Example 1 ---------------------------------------------------*/
      std::cout << "Running Example 1: Basic 1D Uncontrolled" << std::endl;
      uq_solver = std::make_shared<UQSolver1D>(Basic(nx,ns));
      fixed_control_ev     = true;
      cdf_control_strategy = "none";

      /** Compute CDFs */
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);
      break;
    case 2:
      /** ------ Example 2 ---------------------------------------------------*/
      std::cout << "Running Example 2: Different Speeds" << std::endl;
      uq_solver = std::make_shared<UQSolver1D>(SpeedTest(nx,ns));
      fixed_control_ev     = true;
      cdf_control_strategy = "none";

      /** Compute CDFs */
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);
      break;
    case 3:
      /** ------ Example 3 ---------------------------------------------------*/
      std::cout << "Running Example 3: 2D Uncontrolled" << std::endl;
      uq_solver = std::make_shared<UQSolver2D>(TwoDimensional(nx,ns));
      fixed_control_ev     = true;
      cdf_control_strategy = "none";

      /** Compute CDFs */
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);
      break;
    case 4:
      /** ------ Example 4 ---------------------------------------------------*/
      std::cout << "Running Example 4: CDF Bounds" << std::endl;
      uq_solver = std::make_shared<UQSolver1D>(CDFBounds(nx,ns));
      fixed_control_ev     = false;
      cdf_control_strategy = "none";

      /** Compute CDF bounds */
      num_bound = true;
      uq_solver->computeCDFBounds(use_min_cost, cdf_control_strategy, num_bound);

      /** Compute CDFs with asymmetrical lambdas */
      for (int i = 1; i < 5; i++) {
        for (int j = 4; j > 0; j--) {
          uq_solver = std::make_shared<UQSolver1D>(AsymmetricLambda(nx,ns,i,j));
          uq_solver->computeCDF(use_min_cost, cdf_control_strategy);
        }
      }
      break;
    case 5:
      /** ------ Example 5 ---------------------------------------------------*/
      std::cout << "Running Example 5: 1D Controlled" << std::endl;
      fixed_control_ev = false;

      /** CDF Optimization: Compute EV and then CDFs */
      cdf_control_strategy = "CDF";
      uq_solver = std::make_shared<UQSolver1D>(Controlled(nx,ns,cdf_control_strategy));
      uq_solver->computeExpectedValue(fixed_control_ev);
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);

      /** CDF Optimization: Monte Carlo simulations */
      init_location = {0.4};
      init_mode = 0;
      deadlines = {"0.38", "0.84"};
      for (unsigned long i = 0; i < deadlines.size(); i++) {
        mc_control_strategy = "deadline_" + deadlines[i];
        uq_solver->monteCarlo(init_location, init_mode, n_trials, mc_control_strategy);
      }
      uq_solver->monteCarlo(init_location, init_mode, n_trials, "EV");

      /** EV Optimization: Compute EV and then CDFs */
      cdf_control_strategy = "EV";
      uq_solver = std::make_shared<UQSolver1D>(Controlled(nx,ns,cdf_control_strategy));
      uq_solver->computeExpectedValue(fixed_control_ev);
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);

      /** EV Optimization: Monte Carlo simulations */
      mc_control_strategy = "EV";
      uq_solver->monteCarlo(init_location, init_mode, n_trials, mc_control_strategy);
      break;
    case 6:
      /** ------ Example 6 ---------------------------------------------------*/
      std::cout << "Running Example 6: 2D Controlled" << std::endl;
      fixed_control_ev = false;

      /** CDF Optimization: Compute EV and then CDFs */
      cdf_control_strategy = "CDF";
      uq_solver = std::make_shared<UQSolver2D>(TwoDimensionalControlled(nx,ns,cdf_control_strategy));
      uq_solver->computeExpectedValue(fixed_control_ev);
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);

      /** CDF Optimization: Monte Carlo simulations */
      init_location = {0.4, 0.3};
      init_mode = 0;
      deadlines = {"0.28", "0.33", "0.4"};
      for (unsigned long i = 0; i < deadlines.size(); i++) {
        mc_control_strategy = "deadline_" + deadlines[i];
        uq_solver->monteCarlo(init_location, init_mode, n_trials, mc_control_strategy);
      }
      uq_solver->monteCarlo(init_location, init_mode, n_trials, "EV");

      /** EV Optimization: Compute EV and then CDFs */
      cdf_control_strategy = "EV";
      uq_solver = std::make_shared<UQSolver2D>(TwoDimensionalControlled(nx,ns,cdf_control_strategy));
      uq_solver->computeExpectedValue(fixed_control_ev);
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);

      /** EV Optimization: Monte Carlo simulations */
      mc_control_strategy = "EV";
      uq_solver->monteCarlo(init_location, init_mode, n_trials, mc_control_strategy);
      break;
    case 7:
      /** ------ Example 7 ---------------------------------------------------*/
      std::cout << "Running Appendix: Fishing and Harvesting" << std::endl;
      uq_solver = std::make_shared<UQSolver1D>(Fishing(nx,ns));
      fixed_control_ev     = true;
      cdf_control_strategy = "none";
      
      /** Compute CDFs */
      uq_solver->computeCDF(use_min_cost, cdf_control_strategy);
      break;
  }

  /** Stop timer */
  auto t2 = std::chrono::high_resolution_clock::now();
  auto time_elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
  std::cout << "Complete, total execution time was " << time_elapsed << " milliseconds." << std::endl;
}