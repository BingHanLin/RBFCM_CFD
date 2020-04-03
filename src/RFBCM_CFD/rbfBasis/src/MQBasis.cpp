#include "MQBasis.hpp"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

// for nuemann, any better idea?
// void MQBasis::setOperatorStatus(operatorType OOperatorStatus,
//                                   const std::vector<double>& NNormVec){
// inputOperatorType = OOperatorStatus;
// NormVec.resize( NNormVec.size() );
// copy( NNormVec, NormVec );
// };

double MQBasis::getBasisValue(const vec3d<double>& nodeI,
                              const vec3d<double>& nodeJ,
                              const operatorType inputOperatorType)
{
    double rs = (nodeI - nodeJ).squaredNorm();
    vec3d<double> rr = nodeI - nodeJ;

    double temp;
    if (inputOperatorType == operatorType::IdentityOperation)
    {
        temp = std::sqrt(rs + cc_ * cc_);
    }
    else if (inputOperatorType == operatorType::Laplace)
    {
        double temp1;
        temp1 = std::sqrt(rs + cc_ * cc_) * (rs + cc_ * cc_);
        temp = (rs + 2 * cc_ * cc_) / temp1;
    }
    else if (inputOperatorType == operatorType::Partial_D1)
    {
        double temp1;
        temp1 = std::sqrt(rs + cc_ * cc_);
        temp = rr(0) / temp1;
    }
    else if (inputOperatorType == operatorType::Partial_D2)
    {
        double temp1;
        temp1 = std::sqrt(rs + cc_ * cc_);
        temp = rr(1) / temp1;
    }
    else if (inputOperatorType == operatorType::Partial_D3)
    {
        double temp1;
        temp1 = std::sqrt(rs + cc_ * cc_);
        temp = rr(2) / temp1;
    }
    // else if (inputOperatorType == operatorType::NEUMANN_OPERATOR)
    // {
    //     double temp1;
    //     temp1 = std::sqrt(rs + cc_ * cc_);
    //     temp = 0.0;

    //     for (int dim = 0; dim < dim; dim++)
    //     {
    //         temp += NormVec[dim] * rr[dim] / temp1;
    //     }
    // }
    else
    {
        std::cout << "Operator is not defined!" << std::endl;
        assert(false);
    }

    return temp;
}

Eigen::VectorXd MQBasis::collectOnNodes(
    const std::vector<vec3d<double>>& nodesCloud,
    const operatorType inputOperatorType)
{
    const int neighborNum = nodesCloud.size();
    Eigen::MatrixXd phi(neighborNum, neighborNum);

    for (int i = 0; i < neighborNum; i++)
    {
        for (int j = 0; j < neighborNum; j++)
        {
            vec3d<double> nodeI = nodesCloud[i];
            vec3d<double> nodeJ = nodesCloud[j];

            phi(j, i) /* transposed */ =
                getBasisValue(nodeI, nodeJ, operatorType::IdentityOperation);
        }
    }

    Eigen::VectorXd phiL(neighborNum);
    for (int j = 0; j < neighborNum; j++)
    {
        vec3d<double> nodeI = nodesCloud[0];
        vec3d<double> nodeJ = nodesCloud[j];

        phiL(j) = getBasisValue(nodeI, nodeJ, inputOperatorType);
    }

    Eigen::VectorXd x(neighborNum);
    x = phi.ldlt().solve(phiL);

    return x;
}
