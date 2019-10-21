#include "MQBasis2D.hpp"
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>

// for nuemann, any better idea?
// void MQBasis2D::setOperatorStatus(operatorType OOperatorStatus,
//                                   const std::vector<double>& NNormVec){
// inputOperatorType = OOperatorStatus;
// NormVec.resize( NNormVec.size() );
// copy( NNormVec, NormVec );
// };

double MQBasis2D::getBasisValue(const std::vector<double>& nodeI,
                                const std::vector<double>& nodeJ,
                                const operatorType inputOperatorType)
{
    assert(nodeI.size() == nodeJ.size());

    int dim = nodeI.size();

    std::vector<double> rr(dim);
    double rs = 0;
    for (int i = 0; i < dim; i++)
    {
        rr[i] = nodeI[i] - nodeJ[i];
        rs += rr[i] * rr[i];
    }

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
        temp = rr[0] / temp1;
    }
    else if (inputOperatorType == operatorType::Partial_D2)
    {
        double temp1;
        temp1 = std::sqrt(rs + cc_ * cc_);
        temp = rr[1] / temp1;
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

Eigen::VectorXd MQBasis2D::collectOnNodes(
    const std::vector<std::vector<double>>& nodesCloud,
    const operatorType inputOperatorType)
{
    int neighborNum = nodesCloud.size();

    std::vector<double> nodeI(neighborNum);
    std::vector<double> nodeJ(neighborNum);

    Eigen::MatrixXd phi(neighborNum, neighborNum);
    for (int i = 0; i < neighborNum; i++)
    {
        for (int j = 0; j < neighborNum; j++)
        {
            nodeI = nodesCloud[i];
            nodeJ = nodesCloud[j];

            phi(j, i) /* transposed */ =
                getBasisValue(nodeI, nodeJ, inputOperatorType);
        }
    }

    Eigen::VectorXd phiL(neighborNum);
    for (int j = 0; j < neighborNum; j++)
    {
        nodeI = nodesCloud[0];
        nodeJ = nodesCloud[j];

        phiL(j) = getBasisValue(nodeI, nodeJ, operatorType::IdentityOperation);
    }

    Eigen::VectorXd x(neighborNum);
    x = phi.colPivHouseholderQr().solve(phiL);

    return x;
}
