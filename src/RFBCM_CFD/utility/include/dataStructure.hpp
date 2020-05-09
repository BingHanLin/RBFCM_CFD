#ifndef DATASTRUCTURE_HPP
#define DATASTRUCTURE_HPP

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <Eigen/Sparse>

template <typename T>
using vec3d = Eigen::Matrix<T, 3, 1>;

struct nodesCloud
{
    nodesCloud(const int size)
    {
        id.resize(size);
        nodes.resize(size);
    }
    std::vector<int> id;
    std::vector<vec3d<double>> nodes;
};

#endif
