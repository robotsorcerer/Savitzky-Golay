/*  
*   MatrixXi vander(const int F, double A[25])
*     -Computes the vandermonde matrix and flips it column-wise from left to right
*
*   MatrixXd B = MatrixXd sgdiff(int k, double F) 
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
#include <Eigen/Dense>
#include <Eigen/QR>

using namespace Eigen;
using namespace std;

//Global Variables
int F = 5;      //Frame Size
int k = 3;      //Example Polynomial Order

// Function Prototypes
MatrixXi vander(const int F, double A[25]);
MatrixXd sgdiff(int k, double Fd);
void savgolfilt(VectorXf x, int k, int F, MatrixXd DIM);

double A[25];   //Default vandermonde matrix

class RMatrix
{
  public:
    MatrixXd R;

};

MatrixXi vander(const int F, double A[25])
{
  VectorXi v = VectorXi::LinSpaced(F,(-(F-1)/2),((F-1)/2)).transpose().eval();
  //cout << "v: " << v << endl;

  int j;
  int i;
  int k;

  for (j = 0; j < F; j++) {
    i = j * F;
    for (k = 0; k < F; k++) {
      A[i + k] = v[k];
    }
  }

  for (j = 0; j < F; j++) {
    A[20 + j] = 1.0;
  }

  for (j = 0; j < F; j++) {
    for (k = 0; k < 4; k++) {
      A[(j - (k + 1) * F) + 20] *= A[(j - k * F) + 20];
    }
  }

   MatrixXi A_fliplr(F,F);                      //flip the resulting matrix from left to right
   A_fliplr << A[20], A[15], A[10], A[5], A[0], 
                A[21], A[16], A[11], A[6], A[1],
                A[22], A[17], A[12], A[7], A[2],
                A[23], A[18], A[13], A[8], A[3],
                A[24], A[19], A[14], A[9], A[4];

  return A_fliplr;
}

/*Compute S-Golay Matrix of differentiators*/

MatrixXd sgdiff(int k, double Fd)
{
  //We set the weighting matrix to an identity matrix if no weighting matrix is supplied
  MatrixXd W = MatrixXd::Identity(Fd, Fd);      
//  cout << "\nWeighting Matrix: \n" << W << endl;

  //Compute Projection Matrix B
  MatrixXi s = vander(F, A);   

  //Retrieve vandermonde block from projection matrix
  MatrixXi S = s.block(0, 0, s.rows(), (k+1) ) ; 

//  cout << "\nVandermonde Matrix: \n" << S << endl;

  //Compute sqrt(W)*S
  MatrixXd Sd = S.cast<double> ();    //cast S to double
  MatrixXd inter = W * Sd;              //W is assumed to be identity. Change this if you have reasons to.
//  cout << "\nIntermediate matrix: \n" << inter << endl;

  //Compute the QR Decomposition
  HouseholderQR<MatrixXd> qr(inter);
  qr.compute(inter);

  FullPivLU<MatrixXd>lu_decomp(inter);      //retrieve rank of matrix
  
  int Rank = lu_decomp.rank() ;

//  cout <<"\nInter's Rank is: " << Rank << endl;

  //MatrixXd Rnew = RMatrix.R;
/*
  if (Rank == inter.rows() + 1 || Rank == inter.cols() + 1 )
  {
    MatrixXd R = qr.matrixR().template triangularView<Upper>();         //retrieve the R - Matrix
  
    cout << "\nR Matrix: \n" <<  R << endl;   

    //Compute Matrix of Differentiators
    MatrixXd G = Sd * R.inverse();

    cout << "\n Matrix G: \n" << G << endl;
  }

  else    */      
    //For rank deficient matrices        
    MatrixXd Q = qr.householderQ();
    MatrixXd R = qr.matrixQR().topLeftCorner(Rank, Rank).template triangularView<Upper>();
  //  cout << "\n Q is: " << Q << endl;
  //  cout <<"\n R with rank deficiency: \n" << R << endl; 

    //Compute Matrix of Differentiators
    MatrixXd Rinv = R.inverse();
    MatrixXd RinvT = Rinv.transpose();
    MatrixXd G = Sd * Rinv * RinvT;           //G is the matrix of differentiation filters
/*
    cout << "\n R inverse: \n" << Rinv << endl;
    cout << "\n R inv transpose: \n" << RinvT << endl;
    cout << "\n Matrix G: \n" << G << endl;
    */
 
    MatrixXd SdT = Sd.transpose().eval();

    //cout <<"\nInter Transposed: \n" << SdT << endl;

    //Find the matrix of differentiators
    MatrixXd B = G * SdT * W;

    //cout << "\n Matrix B: \n" << B << endl;

    return B;
}

void savgolfilt(VectorXf x, int k, int F)
{  
  Matrix4d DIM = Matrix4d::Zero();        //initialize DIM as a matrix of zeros if it is not supplied
  cout <<"\nDIM: \n" << DIM <<endl;
  //Reshape depth values by working along the first non-singleton dimension

}


int main ()
{
  MatrixXi s = vander(F, A);        //Compute vandermonde matrix

  cout << "\n Vandermonde Matrix: \n" << s  << endl;

  double Fd = (double) F;        //sets the frame size for the savgol differentiation coefficients. This must be odd

  MatrixXd B = sgdiff(k, Fd);

  cout << "\nMatrix of Differentiators: \n" << B << endl;

  VectorXf x = VectorXf::LinSpaced(80, 900, 980);

  savgolfilt(x, k, F);

  return 0;
}

/* Compile:
cd ../; rm build -rf; mkdir build; cd build; cmake ../; make; ./savgol
*/