/*  
*   MatrixXi vander(const int F)
*     -Computes the vandermonde matrix (the polynomial basis vectors) and flips it column-wise from left to right
*
*   MatrixXf B = MatrixXf sgdiff(int k, double F) 
*      - designs a Savitzky-Golay (polynomial) FIR smoothing
*   filter B.  The polynomial order, k, must be less than the frame size of the convolution coefficients,
*   F. F must be odd. 
*
*   Author: Olalekan Ogunmolu  
*           August 12, 2015
*   
    Reference: INTRODUCTION TO SIGNAL PROCESSING [Chapter 8; Section 8.3.5]
                Sophocles J. Orfanidis, Prentice Hall, 2010
*/

// Include Files
#include "savgol.h"
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/QR>
#include <cmath>

using namespace Eigen;
using namespace std;

//Global Variables
int F = 5;      //Frame Size
int k = 3;      //Example Polynomial Order
double Fd = (double) F;        //sets the frame size for the savgol differentiation coefficients. This must be odd

// Function Prototypes
//MatrixXi vander(const int F, double A[25]);
MatrixXi vander(const int F);
MatrixXf sgdiff(int k, double Fd);
void savgolfilt(VectorXf x, VectorXf x_on, int k, int F, MatrixXf DIM);

/*Compute the polynomial basis vectors s_0, s_1, s_2 ... s_n using the vandermonde matrix.
This is super-hacky!  
*/
MatrixXi vander(const int F)
{
  VectorXi v = VectorXi::LinSpaced(F,(-(F-1)/2),((F-1)/2)).transpose().eval();

  MatrixXi A(F, F+1);     //We basically compute an F X F+1 matrix;

  for(int i = 0; i < F; ++ i)
  {
    for(int j=1; j < F+1; ++j)
    {
      A(i,j) = pow(v(i), (j-1) ); 
    }
  }

  A = A.block(0, 1, F, F );   //and retrieve the right F X F, excluding the first column block to find the vandermonde matrix.

  return A;
}

/*Compute the S-Golay Matrix of differentiators*/
MatrixXf sgdiff(int k, double Fd)
{
  //We set the weighting matrix to an identity matrix if no weighting matrix is supplied
  MatrixXf W = MatrixXf::Identity(Fd, Fd);      

  //Compute Projection Matrix B
  MatrixXi s = vander(F);   

  //Retrieve vandermonde block from projection matrix
  MatrixXi S = s.block(0, 0, s.rows(), (k+1) ) ; 

  //Compute sqrt(W)*S
  MatrixXf Sd = S.cast<float> ();    //cast S to float
  MatrixXf inter = W * Sd;              //W is assumed to be identity. Change this if you have reasons to.

  //Compute the QR Decomposition
  HouseholderQR<MatrixXf> qr(inter);
  qr.compute(inter);

  FullPivLU<MatrixXf>lu_decomp(inter);      //retrieve rank of matrix
  
  int Rank = lu_decomp.rank() ;
   
  //For rank deficient matrices. The S matrix block will always be rank deficient.        
  MatrixXf Q = qr.householderQ();
  MatrixXf R = qr.matrixQR().topLeftCorner(Rank, Rank).template triangularView<Upper>();

  //Compute Matrix of Differentiators
  MatrixXf Rinv = R.inverse();
  MatrixXf RinvT = Rinv.transpose();

  MatrixXf G = Sd * Rinv * RinvT;           /*G = S(S'S)^(-1)   -- eqn 8.3.90 (matrix of differentiation filters)*/
 
  MatrixXf SdT = Sd.transpose().eval();

  MatrixXf B = G * SdT * W;   //SG-Smoothing filters of length F and polynomial order k

  return B;
}

void savgolfilt(VectorXf x, VectorXf x_on, int k, int F)
{  
  Matrix4f DIM = Matrix4f::Zero();        //initialize DIM as a matrix of zeros if it is not supplied
  int siz = x.size();       //Reshape depth values by working along the first non-singleton dimension

  //Find leading singleton dimensions
  
  //Pre-allocate output vector
  VectorXf y(siz);

  MatrixXf B = sgdiff(k, Fd);       //retrieve matrix B

  /*Transient On*/
  int id_size = (F+1)/2 - 1;
  //VectorXf y_on = VectorXf::LinSpaced(id_size, 1, (F+1)/2-1) ;    //preallocate y
  MatrixXf Bbutt = B.bottomLeftCorner((F-1)/2, B.cols());

  int n = Bbutt.rows();
  //flip up and down Bbutt
  MatrixXf Bbuttflipped(Bbutt.rows(), Bbutt.cols());
 
    for(int j = n - 1; j >= 0;)
    { 
      for(int i = 0; i < n ; ++i)
      {        
        Bbuttflipped.row(i) = Bbutt.row(j);
        j--;
      }
    }
  //flip x_on up and down as above
  VectorXf x_onflipped(x_on.rows(), x_on.cols());  //pre-allocate
  VectorXf x_onflippedT = x_onflipped.transpose().eval();     //turn x_on to column vector

  int m = x_on.size();                          //retrieve total # coefficients

    for(int j = m -1; j >=0;)
    {
      for(int i = 0; i < m; ++i)
      {
        x_onflippedT.row(i) = x_on.row(j);
        j--;
      }
    }
  
  //Now compute the transient on
  VectorXf y_on = Bbuttflipped * x_onflipped;

 /*Compute the steady state output*/
  size_t idzeroth = floor(B.cols()/2);
  VectorXf Bzeroth = B.col(idzeroth);
  VectorXf Bzerothf = Bzeroth.cast<float>();
  y(0) = Bzerothf.transpose().eval() * x;

  /*Some error checking criteria*/  
  cout << "\nB: \n" << B <<
          "\n\nMiddle Column of B: \n" << Bzeroth <<
          "\nBzeroth dims: " << Bzeroth.rows() <<" X " << Bzeroth.cols() << 
          "\n\nx dims: " << x.rows() <<" X " << x.cols() << 
          "\n\nBbutt: \n" << Bbutt << 
          "\n\nBbuttflipped: \n" << Bbuttflipped << 
          "\n\nx_onflipped: " << x_on <<
          "\n\nx_onflippedT: " << x_onflippedT <<          
          "\ny_on: " << y_on.transpose().eval() <<  
          "\n\nsize of x_onflipped: " << x_onflipped.rows() << ", " << x_onflipped.cols() << endl; 
}


int main ()
{
  MatrixXi s = vander(F);        //Compute vandermonde matrix

  cout << "\n Vandermonde Matrix: \n" << s  << endl;

  MatrixXf B = sgdiff(k, Fd);

  VectorXf x_on = VectorXf::LinSpaced(F, 1, F);     //collect the first five values into a matrix

  //To express as a real filtering operation, we shift x around the nth time instant
  VectorXf x = VectorXf::LinSpaced(5, 900, 980);

  savgolfilt(x, x_on, k, F);

  return 0;
}

/* Compile:
cd ../; rm build -rf; mkdir build; cd build; cmake ../; make; ./savgol
*/