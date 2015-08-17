#ifndef __SAVGOL_H__
#define __SAVGOL_H__

// Include Files
#include <stddef.h>
#include <stdlib.h>
#include "rtwtypes.h"
#include "vander_types.h"
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/Core>
#include <Eigen/QR>

using namespace Eigen;

// Function Declarations
extern MatrixXi vander(const int F);

#endif
