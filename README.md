###INTRODUCTION
 
This code nicely computes the Vandermonde matrix, Savitzky-Golay differentiation filters and smoothing coefficients for a to a noisy, sequantial signal. It is a textbook implementation of the Savitzky-Golay Filter as outlined in Chapter 8 of Sophocles J. Orfanidis's _Intro. to Signal Processing Book_.

####DEPENDENCIES

In order to be able to compile this file, you would need to install the [Eigen3 Library](http://eigen.tuxfamily.org/index.php?title=Main_Page) to do the linear algebra of the the matrices, vectors and related algorithms I used. You can install it by downloading the 3.2.5 library which I used from [here](http://bitbucket.org/eigen/eigen/get/3.2.5.tar.gz) and following the `README` instructions after unpacking the tarball to install.

####USAGE

*  `MatrixXi vander(const int **F**)`
    	-computes the [vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) and the polynomial of basis vectors and flips it column-wise from left to right

*   `MatrixXf B = MatrixXf sgdiff(int **k**, double **F**)`	
		- designs a Savitzky-Golay FIR smoothing filter B with polynomial order _**k**_ and frame size _**F**_ of the convolution coefficients.  The polynomial order, _**k**_, must be less than the frame size **F** and **F** must be odd. 

####Compilation

There is a `CMakeLists.txt` file in the project root folder. From the project root directory:

1.	Create a build directory: `mkdir build && cd build`
2. 	Compile the `cpp` code: 	`cmake ../`
3.	Build your executable: `make`
4. 	Run the executable:	`./savgol`

If you have issues running the files, please use the issues tab to open a bug or issues. I will generally respond within a 24-hour period.
       
 >Reference: **INTRODUCTION TO SIGNAL PROCESSING** [Chapter 8; Section 8.3.5]
                Sophocles J. Orfanidis, Prentice Hall, 2010