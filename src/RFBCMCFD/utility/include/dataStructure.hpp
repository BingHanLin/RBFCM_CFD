#ifndef DATASTRUCTURE_HPP
#define DATASTRUCTURE_HPP

#define EIGEN_INITIALIZE_MATRICES_BY_ZERO

#include <Eigen/Sparse>

template <typename T>
using vec3d = Eigen::Matrix<T, 3, 1>;

struct nodesCloud
{
    size_t id0_;
    std::vector<size_t> ids_;
    size_t size_;

    nodesCloud() : id0_(), ids_()
    {
        size_ = ids_.size();
    }
    nodesCloud(const size_t& id0, const std::vector<size_t>& ids)
        : id0_(id0), ids_(ids)
    {
        size_ = ids_.size();
    }
};

#endif
