/*  
*  Copyright August 2015
*  Author: Olalekan P. Ogunmolu
*
* Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*      http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* 
* See the License for the specific language governing permissions and
* limitations under the License.
* 
*/

// Include Files
#include "savgol.h"
#include <iostream>
#include <cmath>

using namespace Eigen;
using namespace std;

//Global Variables
int F;      //Frame Size
int k;      //Example Polynomial Order
double Fd ;

// Function Prototypes
MatrixXi vander(const int F);
MatrixXf sgdiff(int k, double Fd);
RowVectorXf savgolfilt(VectorXf x, VectorXf x_on, int k, int F, MatrixXf DIM);

/*Compute the polynomial basis vectors s_0, s_1, s_2 ... s_n using the vandermonde matrix.*/
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

  A = A.block(0, 1, F, F );   //and retrieve the right F X F matrix block, excluding the first column block to find the vandermonde matrix.

  return A;
}

/*Compute the S-Golay Matrix of differentiators*/
MatrixXf sgdiff(int k, double Fd)
{
  //We set the weighting matrix to an identity matrix if no weighting matrix is supplied
  MatrixXf W = MatrixXf::Identity(Fd, Fd);      

  //Compute Projection Matrix B
  MatrixXi s = vander(F);   

  //Retrieve the rank deficient matrix from the projection matrix
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

RowVectorXf savgolfilt(VectorXf x, VectorXf x_on, int k, int F)
{  
  Matrix4f DIM = Matrix4f::Zero();        //initialize DIM as a matrix of zeros if it is not supplied
  int siz = x.size();       //Reshape depth values by working along the first non-singleton dimension

  //Find leading singleton dimensions
  
  MatrixXf B = sgdiff(k, Fd);       //retrieve matrix B

  /*Transient On*/
  int id_size = (F+1)/2 - 1;
  MatrixXf Bbutt = B.bottomLeftCorner((F-1)/2, B.cols());

  int n = Bbutt.rows();
  //flip Bbutt from top all the way down 
  MatrixXf Bbuttflipped(n, Bbutt.cols());
 
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
  x_onflipped.transpose().eval();     

  int m = x_on.size();                          //retrieve total # coefficients

    for(int j = m -1; j >=0;)
    {
      for(int i = 0; i < m; ++i)
      {
        x_onflipped.row(i) = x_on.row(j);
        j--;
      }
    }
  
  VectorXf y_on = Bbuttflipped * x_onflipped;  //Now compute the transient on

 /*Compute the steady state output*/
  size_t idzeroth = floor(B.cols()/2);
  VectorXf Bzeroth = B.col(idzeroth);
  VectorXf Bzerothf = Bzeroth.cast<float>();

  VectorXf y_ss = Bzerothf.transpose().eval() * x;     //This is the steady-state smoothed value

  /*Compute the transient off for non-sequential data*/
  MatrixXf Boff = B.topLeftCorner((F-1)/2, B.cols());

  int p = Boff.rows();                        //flip Boff along the horizontal axis

  MatrixXf Boff_flipped(p, Boff.cols());
    
  for(int j = p - 1; j >= 0;)
  { 
    for(int i = 0; i < p ; ++i)
    {        
      Boff_flipped.row(i) = Boff.row(j);
      j--;
    }
  }

/*x_off will be the last (F-1) x values. Note, if you are smoothing in real time, you need to find 
  a way to let your compiler pick the last F-length samples from your data in order to compute your x_off. 
  You could have the program wait for x_milliseconds before you pick 
  the transient off, for example*/
  VectorXf x_off = VectorXf::LinSpaced(F, x(0), x(F-1)).transpose();  
  VectorXf x_offflipped(x_off.rows(), x_off.cols());      //pre-allocate    
  //flip x_off along the horizontal axis
    int q = x_off.size();                          //retrieve total # coefficients

    for(int j = q -1; j >=0;)
    {
      for(int i = 0; i < q; ++i)
      {
        x_offflipped.row(i) = x_on.row(j);
        j--;
      }
    }
  VectorXf y_off = Boff_flipped * x_offflipped;   //This is the transient off

  /*Make Y into the shape of X and retuen the smoothed values!*/
  RowVectorXf y(F);
  y << y_off.transpose().eval(), y_ss, y_on.transpose().eval();

  return y;
}


int main (int argc, char** argv)
{
  float x_min, x_max;
  if(argc>1)
  {    
    for(size_t i = 1; i < argc; ++i )
    {
      F = atoi(argv[1]);
      k = atoi(argv[2]);
      x_min = atoi(argv[3]);
      x_max = atoi(argv[4]);
    }
  }
  else  //use default values
  {
    F = 5; k = 3;
    x_min = 900.0; x_max = 980.0;
  }

  Fd = (double) F;        //sets the frame size for the savgol differentiation coefficients. This must be odd

  MatrixXi s = vander(F);        //Compute vandermonde matrix

  std::cout << "Frame size: " << F << "; \tPolynomial order: " << k << std::endl;
  cout << "\n Vandermonde Matrix: \n" << s  << endl;

  k = atoi(argv[2]) or 3;

  MatrixXf B = sgdiff(k, Fd);

  VectorXf x_on = VectorXf::LinSpaced(F, x_min, x_max);     //collect the first five values into a matrix

  //To express as a real filtering operation, we shift x around the nth time instant
  VectorXf x = VectorXf::LinSpaced(F, x_min, x_max);

  RowVectorXf Filter = savgolfilt(x, x_on, k, F);

  cout <<"\n\nFiltered values in the range \n" << x.transpose().eval() <<"\n are: \n" << Filter << endl;

  return 0;
}

/* Compile:
cd ../; rm build -rf; mkdir build; cd build; cmake ../; make; ./savgol
*/