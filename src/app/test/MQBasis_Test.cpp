
#include "MQBasis.hpp"
#include "catch.hpp"
#include <iostream>

double testFunc(const vec3d<double>& p)
{
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
};

double testFuncD1(const vec3d<double>& p)
{
    return 2 * p[0];
};

TEST_CASE("MQBasis Operator Test", "[RBF]")
{
    SECTION("Derivative 1")
    {
        auto myRBFBasis = std::make_shared<MQBasis>(0.1);

        std::vector<vec3d<double>> nodesD1(10);
        Eigen::VectorXd nodesValue(nodesD1.size());

        for (int nodeID = 0; nodeID < nodesD1.size(); ++nodeID)
        {
            nodesD1[nodeID][0] = nodeID * 0.01;
            nodesValue[nodeID] = testFunc(nodesD1[nodeID]);
        };

        Eigen::SparseMatrix<double> laplaceMatrix(nodesD1.size(),
                                                  nodesD1.size());
        Eigen::SparseMatrix<double> dxMatrix(nodesD1.size(), nodesD1.size());
        Eigen::SparseMatrix<double> dyMatrix(nodesD1.size(), nodesD1.size());
        Eigen::SparseMatrix<double> dzMatrix(nodesD1.size(), nodesD1.size());

        for (int nodeID = 0; nodeID < nodesD1.size(); ++nodeID)
        {
            std::vector<vec3d<double>> neighborNodes(3);
            std::vector<int> neighborID(3);
            if (nodeID == 0)
            {
                neighborNodes[0] = nodesD1[nodeID];
                neighborNodes[1] = nodesD1[nodeID + 1];
                neighborNodes[2] = nodesD1[nodeID + 2];

                neighborID[0] = nodeID;
                neighborID[1] = nodeID + 1;
                neighborID[2] = nodeID + 2;
            }
            else if (nodeID == nodesD1.size() - 1)
            {
                neighborNodes[0] = nodesD1[nodeID - 2];
                neighborNodes[1] = nodesD1[nodeID - 1];
                neighborNodes[2] = nodesD1[nodeID];

                neighborID[0] = nodeID - 2;
                neighborID[1] = nodeID - 1;
                neighborID[2] = nodeID;
            }
            else
            {
                neighborNodes[0] = nodesD1[nodeID - 1];
                neighborNodes[1] = nodesD1[nodeID];
                neighborNodes[2] = nodesD1[nodeID + 1];

                neighborID[0] = nodeID - 1;
                neighborID[1] = nodeID;
                neighborID[2] = nodeID + 1;
            }
            std::cout << neighborID[0] << ", " << neighborID[1] << ", "
                      << neighborID[2] << std::endl;

            Eigen::VectorXd laplaceVector = myRBFBasis->collectOnNodes(
                neighborNodes, rbfOperatorType::LAPLACE);
            Eigen::VectorXd dxVector = myRBFBasis->collectOnNodes(
                neighborNodes, rbfOperatorType::CONSTANT);
            Eigen::VectorXd dyVector = myRBFBasis->collectOnNodes(
                neighborNodes, rbfOperatorType::PARTIAL_D2);
            Eigen::VectorXd dzVector = myRBFBasis->collectOnNodes(
                neighborNodes, rbfOperatorType::PARTIAL_D3);

            for (int i = 0; i < 3; i++)
            {
                laplaceMatrix.insert(nodeID, neighborID[i]) = laplaceVector[i];
                dxMatrix.insert(nodeID, neighborID[i]) = dxVector[i];
                dyMatrix.insert(nodeID, neighborID[i]) = dyVector[i];
                dzMatrix.insert(nodeID, neighborID[i]) = dzVector[i];
            }
        }

        Eigen::VectorXd nodesD1Value = dxMatrix * nodesValue;
        std::cout << dxMatrix << std::endl;

        std::cout << nodesD1Value << std::endl;
        REQUIRE(nodesD1Value[3] == testFunc(nodesD1[3]));
        REQUIRE(nodesD1Value[5] == testFunc(nodesD1[5]));
        REQUIRE(nodesD1Value[8] == testFunc(nodesD1[8]));
    }
}