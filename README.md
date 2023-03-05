# License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------
# Manuscript

The primary purpose in distributing this source code is to enable readers to reproduce the numerical results reported in the manuscript "Quantifying and managing uncertainty in piecewise-deterministic Markov processes" by Elliot Cartee, Antonio Farah, April Nellis, Jacob van Hook, and Alexander Vladimirsky. A pre-print of this article can be found on Arxiv [here](https://arxiv.org/abs/2008.00555). 

--------------------------------------------
# Contributions / Acknowledgements

* TODO: Write this

--------------------------------------------
# Instructions

## Requirements
The C++ code requires an external library:
* [Boost](http://www.boost.org/), which is used for implementation of multidimensional arrays and heaps.

Currently the C++ code is run through a Makefile that assumes this library is installed in the `/usr/local/include/boost` directory.
If this library is installed elsewhere, you will have to modify the [Makefile](https://github.com/eikonal-equation/TimeDependent_SEG/blob/master/Makefile) to make sure the library is properly linked.

The code uses the C++14 standard, and can be compiled using both gcc and icpc.

The plotting scripts are in [Jupyter notebooks](https://jupyter.org), and require the packages [numpy](https://numpy.org/) and [matplotlib](https://matplotlib.org/).

## Running the code

Assuming all necessary libraries/packages are linked, you should be able to compile the code and run the test cases.

To compile, make sure you are in the folder /UQ_PDMP/ and use the command `make` from the command line.

The seven examples in the manuscript can be reproduced using the following commands:
* To compile & run Example 1 (1D, 2 modes: opposite directions with same speed):
` make run TEST='1' `
* To compile & run Example 2 (1D, 2 modes: opposite directions with different speeds):
` make run TEST='2' `
* To compile & run Example 3 (2D, 4 modes: cardinal directions with same speed):
` make run TEST='3' `
* To compile & run Example 4 (1D, 2 modes: bounds on CDFs):
` make run TEST='4' `
* To compile & run Example 5 (1D, 2 modes: controlled):
` make run TEST='5' `
* To compile & run Example 6 (2D, 4 modes: controlled):
` make run TEST='6' `
* To compile & run Example 7 (1D, 3 modes: fish harvesting (appendix)):
` make run TEST='7' `

## Visualizing output

Once you have run the C++ code, you can visualize the output by using the Jupyter notebooks in the directory /UQ_PDMP/plotting/

# Coding conventions

### C++ Naming conventions
The C++ code is written according to the following naming conventions:

* Class names start with a upper case `C` with the rest of the name in Camel case
* Arguments of functions start with a lower case `a`
* Fields start with a lower case `f`
* All functions start with lower case
* Private functions and most inline functions are written in lower case (with underscores)
* All other functions are written in Camel case
* Local variables are in all lower case (with underscores)
* Global constants are in all upper case (with underscores)
* Filenames are in Camel case (except for `main.cpp`)

Otherwise, the code is mostly written to conform to the standards of the [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html)