// Include Files
#include "rt_nonfinite.h"
#include "vander.h"
#include <iostream>
#include <Eigen/LU>
#include <Eigen/Dense>

using namespace std;

// Function Definitions

//
// VANDER Vandermonde matrix.
//    A = VANDER(V) returns the Vandermonde matrix whose columns
//    are powers of the vector V, that is A(i,j) = v(i)^(n-j).
//
//    Class support for input V:
//       float: double, single
// Arguments    : const double v[5]
//                double A[25]
// Return Type  : void
//
void vander(const double v[5], double A[25])
{
  int jtilecol;
  int ibtile;
  int k;

  //    Copyright 1984-2014 The MathWorks, Inc.
  for (jtilecol = 0; jtilecol < 5; jtilecol++) {
    ibtile = jtilecol * 5;
    for (k = 0; k < 5; k++) {
      A[ibtile + k] = v[k];
    }
  }

  for (jtilecol = 0; jtilecol < 5; jtilecol++) {
    A[20 + jtilecol] = 1.0;
  }

  for (jtilecol = 0; jtilecol < 5; jtilecol++) {
    for (k = 0; k < 4; k++) {
      A[(jtilecol - (k + 1) * 5) + 20] *= A[(jtilecol - k * 5) + 20];
    }
  }
  for (int i = 0; i < 25; ++i)
  cout << "A" << "[" << i << "]: " << A[i] << endl;
}

int main ()
{
  double v[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
  double A[25];
  for(int i = 0; i < 5; ++i)
  cout << "v: " << v[i] << endl;
  vander(v, A);
  return 0;
}

/*
cd ../; rm build -rf; mkdir build; cd build; cmake ../
*/