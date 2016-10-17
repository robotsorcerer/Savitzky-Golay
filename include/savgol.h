#ifndef __SAVGOL_H__
#define __SAVGOL_H__

#include <Eigen/Core>

Eigen::MatrixXi vander(const int F);
Eigen::MatrixXf sgdiff(int k, double Fd, int F);
Eigen::RowVectorXf savgolfilt(Eigen::VectorXf const & x, Eigen::VectorXf const & x_on, int k, int F, double Fd);

#endif
