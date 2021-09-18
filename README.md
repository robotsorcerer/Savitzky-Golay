## Savitzky-Golay Filter in C++

Author: [Lekan Ogunmolu](http://ecs.utdallas.edu/~opo140030)
# Table of Contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Compilation](#compilation)
- [Citation](#citation)
- [Issues](#issues)
- [Changelog](#changelog)


### Introduction

Nicely computes the Vandermonde matrix, Savitzky-Golay differentiation filters, and smoothing coefficients for any sequential signal. It is a textbook implementation of the Savitzky-Golay Filter. Initial testing of this code was on a Ubuntu 14.04.02 Trusty OS running Linux 4.4 but will work on any other Linux/Windows/Mac OS machine with little effort.

Below are examples of how the filter smoothes out a noisy depth map data from the kinect time-of-flight sensor:

<img src="/images/Protonect_Uncalibrated.jpg" height="500px" >
<img src="/images/ROS_Calibrated.jpg" height="500px" >
<!-- <img src="/images/Savitzky-Golay_smoothing_filter.jpg" height="500px"> -->

### Dependencies

+  [Eigen3 Library](http://eigen.tuxfamily.org/index.php?title=Main_Page) for the linear algebra of the matrices, vectors and related algorithms. You can download the 3.2.5 library which I used from [here](http://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz) and follow the `README` instructions after unpacking the tarball to install. 

### Usage

* `./savgol`
  
  __Options__

  - -h or --help: print out the help menu.

  - Example: Compute the savitzky-golay filter coefficients with frame size, `F = 5` and polynomial order 3 (these are the default parameters of the filter) for linearly spaced data points between `x_min = 900` and `x_max = 980`. A typical cmd line usage:

    `./savgoal 9 5` where `F = 9` and `k = 5`.

  - To pass in your arbitrary data points between a value `x_min` and `x_max`, run this way in order: ./savgoal `F` `k` `x_min` `x_max`.

  The filtered values are returned to the console. Note that the Frame size should ideally be odd.

### Components
*  `MatrixXi vander(const int F);`
    	
    - Computes the [Vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) and the polynomial of basis vectors; it flips the vector column-wise, left to right.

*   `MatrixXf B = MatrixXf sgdiff(int k, double F)`	
		
	- Designs a Savitzky-Golay FIR smoothing filter, B, with polynomial order _**k**_ and frame size _**F**_ of convolution coefficients.  The polynomial order, _**k**_, must be less than the frame size _**F**_, and _**F**_ must be odd. 

*   `savgolfilt(x, x_on, k, F)`
	
	- Computes the smoothed values of the signal x, whose tansient on is `x_on` initialized with size F.

*    **Note**
	In calculating the transient off, `x_off` will be the last `(F-1)` `x` values, where `x`'s are the data sequence we want to filter.If you are smoothing data offline, then this code will work seamlessly. Just load your data in the `main()` function where, for an example, I have used linearly spaced values between `900` and `980` at a frame `5` size for my steady state values. 
	
	Note, if you are smoothing data in real time, you need to find a way to let your compiler pick the last F-length samples from your data in order to compute your transient off, i.e., x_off. You could have the program wait for x_milliseconds after stopping your code before you pick the transient off, for example.

### Compilation

There is a `CMakeLists.txt` file in the project root folder. From the project root directory:

1.	Create a build directory: `mkdir build && cd build`
2. 	Compile the `cpp` code: 	`cmake ../`
3.	Build your executable: `make`
4. 	Run the executable:	`./savgol`


### Citation

If you have used `Savitzky-Golay` in your work, please cite it.

```tex
@misc{Savitzky-Golay,
  author = {Ogunmolu, Olalekan},
  title = {{Savitzky-Golay Filter in C++}},
  year = {2015},
  howpublished = {\url{https://github.com/lakehanne/Savitzky-Golay}},
  note = {Accessed August 15, 2015}
}
```
#### Issues

If you have issues running the files, please use the issues tab to open a bug. I will generally respond within a 24-hour period.

### Changelog
*	Added citation to README (August 14, 2015)
* Added examples to `int main()` function (August 15, 2015)
* Modified frame size and polynomial order to be reconfigurable at run time (July 1, 2016)
       
### TODO
Add a plotter to plot the filtered values on a gtk chart?

### Reference

**INTRODUCTION TO SIGNAL PROCESSING** 

  Sophocles J. Orfanidis, Prentice Hall, 2010

  *Chapter 8; Section 8.3.5*
