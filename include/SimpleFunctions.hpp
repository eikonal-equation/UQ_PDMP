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
* File: SimpleFunctions.hpp
*
* Authors: Elliot Cartee, April Nellis, Jacob Van Hook
*
* Description: This file contains inline definitions of many small/simple
* functions used for speeds or running costs.
*
* ==============================================================================
*/
#ifndef SIMPLEFUNCTIONS_HPP
#define SIMPLEFUNCTIONS_HPP

/** ------ Libraries ---------------------------------------------------------*/
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConstants.hpp"

/** Velocity functions */
namespace dynamics {

/** ------ One-dimensional velocity functions --------------------------------*/
inline double constant_right(double x, double a) {
  return 1.0;
}

inline double constant_left(double x, double a) {
  return -1.0;
}

inline double slow_right(double x, double a) {
  return 0.5;
}

inline double control_right(double x, double a) {
  return 0.5 + a;
}

inline double control_left(double x, double a) {
  return -0.5 + a;
}

/** Allee-effect fishing population dynamics
  * Current fish population represented by x */
inline double fish_low(double x, double a) {
  return 2 * x * (x - 1) * (1 - (x / 3.8))  -  1.057 * x;
}

inline double fish_mid(double x, double a) {
  return 2 * x * (x - 1) * (1 - (x / 4.0))  -  1.057 * x;
}

inline double fish_high(double x, double a) {
  return 2 * x * (x - 1) * (1 - (x / 4.2))  -  1.057 * x;
}

/** ------ Two-dimensional velocity functions --------------------------------*/
inline double zero(double x, double y, double a) {
  return 0.0;
}

inline double positive(double x, double y, double a) {
  return 1.0;
}

inline double negative(double x, double y, double a) {
  return -1.0;
}

inline double control_2d_right_x(double x, double y, double a) {
  return 0.5 + cos(a);
}

inline double control_2d_right_y(double x, double y, double a) {
  return sin(a);
}

inline double control_2d_left_x(double x, double y, double a) {
  return -0.5 + cos(a);
}

inline double control_2d_left_y(double x, double y, double a) {
  return sin(a);
}

inline double control_2d_up_x(double x, double y, double a) {
  return cos(a);
}

inline double control_2d_up_y(double x, double y, double a) {
  return 0.5 + sin(a);
}

inline double control_2d_down_x(double x, double y, double a) {
  return cos(a);
}

inline double control_2d_down_y(double x, double y, double a) {
  return -0.5 + sin(a);
}

} // namespace speeds

/** ------ Cost functions ----------------------------------------------------*/
namespace costs {

inline double constant2(double x, double a) {
  return 1.0;
}

inline double constant3(double x, double y, double a) {
  return 1.0;
}

/** Portion of fish population being harvested
  * Fish population represented by x */
inline double harvest_cost(double x, double a) {
  return x * 1.057;
}

} // namespace costs

/** Small helper method for reading a command line argument in as an int */
int read_arg(char* argv) {
std::string arg = argv;
  int x = 1;
  try {
    std::size_t pos;
    x = std::stoi(arg, &pos);
    if (pos < arg.size()) {
      std::cerr << "Trailing characters after number: " << arg << '\n';
    }
  } catch (std::invalid_argument const &ex) {
    std::cerr << "Invalid number: " << arg << '\n';
  } catch (std::out_of_range const &ex) {
    std::cerr << "Number out of range: " << arg << '\n';
  }
  return x;
}

#endif