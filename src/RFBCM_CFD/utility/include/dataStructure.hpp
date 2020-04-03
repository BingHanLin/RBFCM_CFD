#ifndef DATASTRUCTURE_HPP
#define DATASTRUCTURE_HPP

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <Eigen/Sparse>

template <typename T>
using vec3d = Eigen::Matrix<T, 3, 1>;

#endif
