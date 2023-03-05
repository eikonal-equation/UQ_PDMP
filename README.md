# License

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <http://www.gnu.org/licenses/>.

--------------------------------------------
# Manuscript

The primary purpose in distributing this source code is to enable readers to reproduce the numerical results reported in the manuscript "Quantifying and managing uncertainty in piecewise-deterministic Markov processes" by Elliot Cartee, Antonio Farah, April Nellis, Jacob van Hook, and Alexander Vladimirsky. A pre-print of this article can be found on Arxiv [here](https://arxiv.org/abs/2008.00555). 

--------------------------------------------
# Contributions / Acknowledgements

* Elliot Cartee was primarily responsible for the overall structure of the final implementation in this repository.
* April Nellis and Jacob van Hook made significant contributions to the computations and plotting of the CDF bounds (Section 3, Example 4) and the fishing example (in the Appendix of our paper).
  During the summer REU program, they were also responsible for an earlier (proof-of concept, 1D) implementation of methods in sections 2 and 4.
* Antonio Farah's main focus was on computing and analyzing the CDF for the number of random mode-switches experienced up to the termination.  
  While interesting in its own right, this topic was omitted from the final manuscript for the sake of brevity & to streamline our presentation.
* Alexander Vladimirsky supervised the overall project, formulated this model of uncertainty quantification, contributed several algorithmic ideas, and suggested the specific examples for benchmarking.
* All five authors have contributed to writing and editing the manuscript.  

The authors would like to thank Tristan Reynoso and Shriya Nagpal for their help in the initial stages of this project during the summer
REU-2018 program at Cornell University.

This work was started in an REU program partially supported by the NSF-RTG award DMS-1645643. 
Elliot Cartee's and Alexander Vladimirsky's work was also supported in part by the NSF award DMS-1738010. 
Alexander Vladimirsky's work was also supported by the Simons Foundation Fellowship and by the NSF award DMS-2111522.

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
