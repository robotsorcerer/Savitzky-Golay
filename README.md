#Savitzky-Golay Filter in C++

##Author: [Olalekan P. Ogunmolu](http://lakehanne.github.io)<<olalekan.ogunmolu@utdallas.edu>>, [SEnsing, Robotics, Vision, Control and Estimation (SERVICE) Lab](http://ecs.utdallas.edu/research/researchlabs/service-lab/), University of Texas at Dallas, Richardson, TX, USA

Copyright August 2015

##Table of Contents
- [CHANGELOG](#changelog)
- [Introduction](#introduction)
- [DEPENDENCIES](#dependencies)
- [USAGE](#usage)
- [COMPILATION AND RUNNING THE PROJECT](#compilation-and-running-the-project)
- [Citation](#citation)
- [ISSUES](#issues)

##CHANGELOG
*	Added citation to README (August 14, 2015)
*   Added examples to `int main()` function (August 15, 2015)

##INTRODUCTION
This code nicely computes the Vandermonde matrix, Savitzky-Golay differentiation filters and smoothing coefficients for any noisy, and sequantial signal. It is a textbook implementation of the Savitzky-Golay Filter. It can be run offline with collected data but needs slight [tweaking](#usage) if it must be run with on-line data in real time. Initial testing of this code was done using a Ubuntu 14.04.02 Trusty OS running Linux 4.4 but will work on any other Linux/Windows/Mac OS machine with little effort, I presume.

##DEPENDENCIES

In order to be able to compile this file, you would need to install the [Eigen3 Library](http://eigen.tuxfamily.org/index.php?title=Main_Page) to do the linear algebra of the matrices, vectors and related algorithms I used. You can install it by downloading the 3.2.5 library which I used from [here](http://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz) and following the `README` instructions after unpacking the tarball to install. I have not tested it with other versions of the library. If you are having problems running this code in other versions of Eigen, please raise an issue using the link on the right of this page or go through the [eigen documentation page](http://eigen.tuxfamily.org/dox/index.html) if you are impatient.

##USAGE

*  `MatrixXi vander(const int F)`
    	
    - computes the [vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) and the polynomial of basis vectors and flips it column-wise from left to right

*   `MatrixXf B = MatrixXf sgdiff(int k, double F)`	
		
	- designs a Savitzky-Golay FIR smoothing filter B with polynomial order _**k**_ and frame size _**F**_ of the convolution coefficients.  The polynomial order, _**k**_, must be less than the frame size _**F**_ and _**F**_ must be odd. 


*	**Note**
	In calculating the transient off, x_off will be the last (F-1) x values, where x's are the data sequence we want to filter.If you are smoothing data offline, then this code will work seamlessly. Just load your data in the `main()` function where, for an example, I have used linearly spaced values between 900 and 980 at a frame 5 size for my steady state values and initiated the transient on with frame sized linearly spaced values between 960 and 980. 
	
	Note, if you are smoothing data in real time, you need to find a way to let your compiler pick the last F-length samples from your data in order to compute your transient off, i.e., x_off. You could have the program wait for x_milliseconds after stopping your code before you pick the transient off, for example.

##COMPILATION AND RUNNING THE PROJECT

There is a `CMakeLists.txt` file in the project root folder. From the project root directory:

1.	Create a build directory: `mkdir build && cd build`
2. 	Compile the `cpp` code: 	`cmake ../`
3.	Build your executable: `make`
4. 	Run the executable:	`./savgol`


##CITATION

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
##ISSUES
If you have issues running the files, please use the issues tab to open a bug. I will generally respond within a 24-hour period.
       
 >Reference: **INTRODUCTION TO SIGNAL PROCESSING** [Chapter 8; Section 8.3.5]
                Sophocles J. Orfanidis, Prentice Hall, 2010