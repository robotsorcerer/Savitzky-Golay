#include "rt_nonfinite.h"
#include "vander.h"
#include "main.h"
#include "vander_terminate.h"
#include "vander_initialize.h"

// Function Declarations
static void argInit_1x5_real_T(double result[5]);
static double argInit_real_T();
static void main_vander();

// Function Definitions

//
// Arguments    : double result[5]
// Return Type  : void
//
static void argInit_1x5_real_T(double result[5])
{
  int b_j1;

  // Loop over the array to initialize each element.
  for (b_j1 = 0; b_j1 < 5; b_j1++) {
    // Set the value of the array element.
    // Change this value to the value that the application requires.
    result[b_j1] = argInit_real_T();
  }
}

//
// Arguments    : void
// Return Type  : double
//
static double argInit_real_T()
{
  return 0.0;
}

//
// Arguments    : void
// Return Type  : void
//
static void main_vander()
{
  double dv0[5];
  double A[25];

  // Initialize function 'vander' input arguments.
  // Initialize function input argument 'v'.
  // Call the entry-point 'vander'.
  argInit_1x5_real_T(dv0);
  vander(dv0, A);
}

//
// Arguments    : int argc
//                const char * const argv[]
// Return Type  : int
//
int main(int, const char * const [])
{
  // Initialize the application.
  // You do not need to do this more than one time.
  vander_initialize();

  // Invoke the entry-point functions.
  // You can call entry-point functions multiple times.
  main_vander();

  // Terminate the application.
  // You do not need to do this more than one time.
  vander_terminate();
  return 0;
}

//g++ main.cpp rtGetInf.cpp rtGetNaN.cpp rt_nonfinite.cpp vander.cpp vander_initialize.cpp vander_terminate.cpp -o main
//
