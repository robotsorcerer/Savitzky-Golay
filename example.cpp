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

#include <Eigen/Core>

#include <iostream>

using namespace Eigen;
using namespace std;

void help()
{
  std::cout << "================================================== \n" 
            << "USAGE:                                             \n"
            << "\n"
            << "./savgol [<frame_size> [<polynomial_order> [<x_low> <x_high>] ] ] \n" 
            << "\n"
            << "       <frame_size>: odd int and ideally greater than\n" 
            << "       <polynomial_order>:  an integer             \n" 
            << "\n"
            << "       <x_low>, <x_high>: min and max limits of    \n"
            << "       linspaced vector to be filtered.            \n"
            << "===================================================\n" 
            << "       Example: ./savgol 9 5   \n\n";;
}


int main (int argc, char** argv)
{
  int F;      //Frame Size
  int k;      //Example Polynomial Order
  double Fd ;
  float x_min, x_max;

  if(argc>1)
  { 
    if(argv[1] == "-h" || "-help")
    {
      help();
      return EXIT_SUCCESS;
    }
    else
    {
      for(auto i = 1; i < argc; ++i )
      {
        F = atoi(argv[1]);
        k = atoi(argv[2]);
        x_min = atoi(argv[3]);
        x_max = atoi(argv[4]);
      }
    }
  }
  else  //use default values
  {
    F = 5; k = 3;
    x_min = 900.0; x_max = 980.0;
  }

  auto s = vander(F);        //Compute vandermonde matrix

  cout << "Frame size: " << F << "; \tPolynomial order: " << k << endl;
  cout << "\n Vandermonde Matrix: \n" << s  << endl;

  k = atoi(argv[2]) or 3;

  auto B = sgdiff(k, F, Fd);

  //To express as a real filtering operation, we shift x around the nth time instant
  auto x = VectorXf::LinSpaced(F, x_min, x_max);

  auto Filter = savgolfilt(x, k, F);

  cout <<"\n\nFiltered values in the range \n" << x.transpose().eval() <<"\n are: \n" << Filter << endl;

  return 0;
}

/* Compile:
cd ../; rm build -rf; mkdir build; cd build; cmake ../; make; ./savgol
*/
