# Alfonsless
This repository contains matlab software that provides an interface for the SOS programming solver [alfonso](https://arxiv.org/abs/1712.01792). The interface looks similar to the [spotless](https://github.com/spot-toolbox/spotless) solver, but the functionality is slightly reduced. The memory usage of this solver is significantly better than SDP based solvers as dimension grows (see the paper linked above for details).

**As of March 7, 2019, the way we are passing problems to the Alfonso solver and generating bases for the cones that the solver uses results in the solver terminating with an inconclusive result on certain large problems.**

## Dependencies
The file `setup.m` in the base directory is referenced by the example problems and provides an outline of which directories in the dependencies need to be added to the path. If there is a more idiomatic way to handle the path variable in Matlab, either make an issue that explains best practices to me, shoot [me](mailto:owhughes@umich.edu) an email, or ideally make a PR. 

* [Chebfun](https://github.com/chebfun/chebfun)
* [Alfonso](https://github.com/dpapp-github/alfonso)
* [partitions](https://www.mathworks.com/matlabcentral/fileexchange/12009-partitions-of-an-integer)
* [Padua2DM](http://www.netlib.org/numeralgo/): Search for Padua2DM to find the download link.
* [spotless](https://github.com/spot-toolbox/spotless) this is included because I use their polynomial utilities. Eventually if this succeeds, spotless 
## Basic usage:
We will illustrate the basic usage of the program by coding up the problem:   
$minimize \int_{[-1, 2]} w(x) dx    
s.t. (u^3 + 3)w_x(u) \geq 0 for u \in [-1, 2]  
      w(u) \geq 1 for u \in [-1, 2]$

where w_x is the derivative of w with respect to x. You can see the full code in `sample_problems/SOS_derivative.m`. 

**Run the program to see an illustration of the major functionality that Alfonsless provides.**

