#ifndef __SAVGOL_H__
#define __SAVGOL_H__

#include <Eigen/Core>

Eigen::MatrixXi vander(const int F);
Eigen::MatrixXf sgdiff(int k, int F, double Fd);
Eigen::RowVectorXf savgolfilt(Eigen::VectorXf const & x, int k, int F);

#endif
