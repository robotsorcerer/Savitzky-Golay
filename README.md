##Savitzky-Golay Filter in C++

Author: [Olalekan P. Ogunmolu](http://lakehanne.github.io)

##Table of Contents
- [Introduction](#introduction)
- [Dependencies](#dependencies)
- [Usage](#usage)
- [Compilation](#compilation)
- [Citation](#citation)
- [Issues](#issues)
- [Changelog](#changelog)


###Introduction
This code nicely computes the Vandermonde matrix, Savitzky-Golay differentiation filters and smoothing coefficients for any noisy, and sequantial signal. It is a textbook implementation of the Savitzky-Golay Filter. Initial testing of this code was done using a Ubuntu 14.04.02 Trusty OS running Linux 4.4 but will work on any other Linux/Windows/Mac OS machine with little effort.

Below are examples of how the filter smoothes out a noisy depth map data from the kinect time-of-flight sensor:

<img src="/include/Protonect_Uncalibrated.jpg" height="500px" >
<img src="/include/ROS_Calibrated.jpg" height="500px" >
<img src="/include/Savitzky-Golay_smoothing_filter.jpg" height="500px">

###Dependencies

In order to be able to compile this file, you would need to install the [Eigen3 Library](http://eigen.tuxfamily.org/index.php?title=Main_Page) to do the linear algebra of the matrices, vectors and related algorithms I used. You can install it by downloading the 3.2.5 library which I used from [here](http://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz) and following the `README` instructions after unpacking the tarball to install. I have not tested it with other versions of the library. If you are having problems running this code in other versions of Eigen, please raise an issue using the link on the right of this page or go through the [eigen documentation page](http://eigen.tuxfamily.org/dox/index.html) if you are impatient.

###Usage

* `./savgol`

  - computes the savitzky-golay filter coefficients with frame size, `F = 5` and polynomial order 3 (these are the default parameters of the filter) for linearly spaced data points between `x_min = 900` and `x_max = 980`.

  - To change the values of the frame size and polynomial order, do

  `./savgoal 9 5` where `F = 9` and `k = 5`.

  - To pass in your arbitrary data points between a value x_min and x_max, pass in the following argum,ents in order: ./savgoal `F` `k` `x_min` `x_max`.

  The filtered values are returned onto the console. Note that the Frame size should ideally be odd

### RESULTS
The savgol filter tries to compute the moving average of the time-series data fed into it. For example, with a frame size of `9` and polynomial order of `5` for numbers linearly spaced between `100` and `1000`, we obtain the following results by running this code:

```bash
Frame size: 9; Polynomial order: 5 

 Vandermonde Matrix: 
     1     -4     16    -64    256  -1024   4096 -16384  65536
     1     -3      9    -27     81   -243    729  -2187 134737
     1     -2      4     -8     16    -32     64   -128      0
     1     -1      1     -1      1     -1      1     -1      1
     1      0      0      0      0      0      0      0      0
     1      1      1      1      1      1      1      1      1
     1      2      4      8     16     32     64    128    256
     1      3      9     27     81    243    729   2187   6561
     1      4     16     64    256   1024   4096  16384  65536

Filtered values in the range 
  100 212.5   325 437.5   550 662.5   775 887.5  1000
 are: 
662.5   775 887.5  1000   550   100 212.5   325 437.5

```

Or for numbers linearly spaced between `100` and `300`, with `F = 7` and `k = 5`, we obtain:

```bash
Frame size: 7; polynomial order: 5 

 Vandermonde Matrix: 
     1     -3      9    -27     81   -243    729
     1     -2      4     -8     16    -32 134865
     1     -1      1     -1      1     -1      0
     1      0      0      0      0      0      0
     1      1      1      1      1      1      1
     1      2      4      8     16     32     64
     1      3      9     27     81    243    729


Filtered values in the range 
    100 133.333 166.667     200 233.333 266.667     300
 are: 
233.333 266.667     300     200     100 133.333 166.667
```

###Components
*  `MatrixXi vander(const int F);`
    	
    - computes the [vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) and the polynomial of basis vectors and flips it column-wise from left to right

*   `MatrixXf B = MatrixXf sgdiff(int k, double F);`	
		
	- designs a Savitzky-Golay FIR smoothing filter B with polynomial order _**k**_ and frame size _**F**_ of the convolution coefficients.  The polynomial order, _**k**_, must be less than the frame size _**F**_ and _**F**_ must be odd. 

*	`savgolfilt(x, x_on, k, F);`
	
	- computes the smoothed values of the signal x, whose tansient on is `x_on` initialized with size F.

*	**Note**
	In calculating the transient off, `x_off` will be the last `(F-1)` `x` values, where `x`'s are the data sequence we want to filter.If you are smoothing data offline, then this code will work seamlessly. Just load your data in the `main()` function where, for an example, I have used linearly spaced values between `900` and `980` at a frame `5` size for my steady state values. 
	
	Note, if you are smoothing data in real time, you need to find a way to let your compiler pick the last F-length samples from your data in order to compute your transient off, i.e., x_off. You could have the program wait for x_milliseconds after stopping your code before you pick the transient off, for example.

###Compilation

There is a `CMakeLists.txt` file in the project root folder. From the project root directory:

1.	Create a build directory: `mkdir build && cd build`
2. 	Compile the `cpp` code: 	`cmake ../`
3.	Build your executable: `make`
4. 	Run the executable:	`./savgol`


###Citation

If you used `Savitzky-Golay` for your work, please cite it.

```tex
@misc{Savitzky-Golay,
  author = {Ogunmolu, Olalekan},
  title = {{Savitzky-Golay Filter in C++}},
  organization = {SEnsing, Robotics, Vision, Control and Estimation Lab},
  address = {University of Texas at Dallas},
  year = {2015},
  howpublished = {\url{https://github.com/lakehanne/Savitzky-Golay}},
  note = {Accessed August 15, 2015}
}
```
####Issues
If you have issues running the files, please use the issues tab to open a bug. I will generally respond within a 24-hour period.

###Changelog
*	Added citation to README (August 14, 2015)
* Added examples to `int main()` function (August 15, 2015)
* Modified frame size and polynomial order to be reconfigurable at run time (July 1, 2016)
       
###TODO
Add a plotter to plot the filtered values on a gtk chart?

###Reference

**INTRODUCTION TO SIGNAL PROCESSING** 

  Sophocles J. Orfanidis, Prentice Hall, 2010

  *Chapter 8; Section 8.3.5*
