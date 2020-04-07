#include "simulationDomain.hpp"
#include "MQBasis.hpp"
#include "rectangle.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iomanip>
#include <iostream>

SimulationDomain::SimulationDomain(std::shared_ptr<MeshData> mesh,
                                   std::shared_ptr<MQBasis> RBFBasis,
                                   nlohmann::json& param)
    : myMesh_(mesh),
      myRBFBasis_(RBFBasis),
      myParams_(param),
      viscous_(0.0),
      density_(0.0),
      tStepSize_(0.0),
      endTime_(0.0),
      crankNicolsonEpsilon_(0.001),
      crankNicolsonMaxIter_(10),
      neighborNum_(0),
      kdTree_()
{
    setUpSimulation();
    showSummary();

    // kdTree_ = KDTreeTsaiAdaptor<std::vector<std::vector<double>>, double, 2>(
    //     myMesh_.getNodes());

    // assembleMatrix();
}

void SimulationDomain::setUpSimulation()
{
    neighborNum_ = myParams_.at("SolverConstrol").at("NeighborNumber");
    tStepSize_ = myParams_.at("SolverConstrol").at("TimeStepSize");
    endTime_ = myParams_.at("SolverConstrol").at("EndTime");

    solverType_ = myParams_.at("SolverConstrol").at("solverType");

    if (solverType_ == solverTypeEnum::NAVIERSTOKES)
    {
        viscous_ = myParams_.at("PhysicsControl").at("Viscosity");
        density_ = myParams_.at("PhysicsControl").at("Density");

        crankNicolsonEpsilon_ = myParams_.at("SolverConstrol")
                                    .at("NavierStokesOptions")
                                    .at("CrankNicolsonEpsilon");

        crankNicolsonMaxIter_ = myParams_.at("SolverConstrol")
                                    .at("NavierStokesOptions")
                                    .at("CrankNicolsonMaxIter");
    }
    else if (solverType_ == solverTypeEnum::POISSON)
    {
    }
}

void SimulationDomain::showSummary()
{
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
    std::cout << std::setfill(' ') << std::setw(50) << "Summary of simulation"
              << std::endl;
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::setfill(' ')
              << std::endl;

    std::cout << "Number of nodes: " << std::setw(8) << myMesh_->numOfNodes()
              << std::endl;
    std::cout << "Time step size: " << std::setw(8) << tStepSize_ << std::endl;
    std::cout << "End time: " << std::setw(8) << endTime_ << std::endl;
    std::cout << "Neighbor number: " << std::setw(8) << neighborNum_
              << std::endl;

    if (solverType_ == solverTypeEnum::NAVIERSTOKES)
    {
        std::cout << "Density:" << std::setw(8) << density_ << std::endl;
        std::cout << "Viscosity:" << std::setw(8) << viscous_ << std::endl;
        std::cout << "CrankNicolson epsilon: " << std::setw(8)
                  << crankNicolsonEpsilon_ << std::endl;
        std::cout << "CrankNicolson max iter: " << std::setw(8)
                  << crankNicolsonMaxIter_ << std::endl;
    }
    else if (solverType_ == solverTypeEnum::POISSON)
    {
    }
}

void SimulationDomain::assembleCoeffMatrix()
{
    if (solverType_ == solverTypeEnum::NAVIERSTOKES)
    {
    }
    else if (solverType_ == solverTypeEnum::POISSON)
    {
    }
}

void SimulationDomain::solveDomain()
{
    if (solverType_ == solverTypeEnum::NAVIERSTOKES)
    {
    }
    else if (solverType_ == solverTypeEnum::POISSON)
    {
    }
}

// void SimulationDomain::assembleMatrix()
// {
//     rhs_ = Eigen::VectorXd::Zero(myMesh_.numAllNodes_);

//     // using size_t for compatibility reasonwith nanoflann.hpp
//     std::vector<std::vector<double>> nodesCloud(neighborNum_);
//     std::vector<size_t> neighbours(neighborNum_);
//     std::vector<double> outDistSqr(neighborNum_);

//     Eigen::VectorXd localVector = Eigen::VectorXd::Zero(neighborNum_);

//     systemVarMatrix_.resize(myMesh_.numAllNodes_, myMesh_.numAllNodes_);
//     systemVarMatrix_.reserve(
//         Eigen::VectorXi::Constant(myMesh_.numAllNodes_, neighborNum_));

//     // go through all interior nodes
//     for (int i = 0; i < myMesh_.numInnNodes_; i++)
//     {
//         localVector = Eigen::VectorXd::Zero(neighborNum_);

//         // use kdtree find indexes of neighbor nodes
//         kdTree_.query(i, neighborNum_, &neighbours[0], &outDistSqr[0]);

//         // store nodes cloud in vector
//         for (int j = 0; j < neighborNum_; j++)
//         {
//             nodesCloud[j] = myMesh_.getNodes()[neighbours[j]];
//         }

//         localVector += myRBFBasis_.collectOnNodes(
//             nodesCloud, RBFBasisType::operatorType::Laplace);

//         for (int j = 0; j < neighborNum_; j++)
//         {
//             systemVarMatrix_.insert(i, neighbours[j]) = localVector(j);
//         }
//     }

//     // go through all boundary nodes
//     for (int i = myMesh_.numInnNodes_; i < myMesh_.numAllNodes_; i++)
//     {
//         localVector = Eigen::VectorXd::Zero(neighborNum_);

//         // use kdtree find indexes of neighbor nodes
//         kdTree_.query(i, neighborNum_, &neighbours[0], &outDistSqr[0]);

//         // store nodes cloud in vector
//         for (int j = 0; j < neighborNum_; j++)
//         {
//             nodesCloud[j] = myMesh_.getNodes()[neighbours[j]];
//         }

//         localVector += myRBFBasis_.collectOnNodes(
//             nodesCloud, RBFBasisType::operatorType::IdentityOperation);

//         for (int j = 0; j < neighborNum_; j++)
//         {
//             systemVarMatrix_.insert(i, neighbours[j]) = localVector(j);
//         }
//     }

//     // go through all boundary nodes,
//     // change this to dirichlet boundary
//     // later
//     for (int i = myMesh_.numInnNodes_; i < myMesh_.numInnNodes_ + 30; i++)
//     {
//         rhs_(i) = 100;
//     }

//     for (int i = myMesh_.numInnNodes_; i < myMesh_.numAllNodes_; i++)
//     {
//         // use kdtree find indexes of neighbor nodes
//         kdTree_.query(i, neighborNum_, &neighbours[0], &outDistSqr[0]);

//         for (int j = 1; j < neighborNum_; j++)
//         {
//             rhs_(neighbours[j]) =
//                 rhs_(neighbours[j]) -
//                 systemVarMatrix_.coeff(i, neighbours[j]) * rhs_(i);
//             systemVarMatrix_.coeffRef(i, neighbours[j]) = 0.0;
//         }
//     }

//     systemVarMatrix_.makeCompressed();

//     Eigen::SparseMatrix<double> systemVarMatrix_adj =
//         systemVarMatrix_.adjoint();

//     std::cout << "\tIS SELFADJOINT: "
//               << (systemVarMatrix_adj.isApprox(systemVarMatrix_) ? "YES\n"
//                                                                  : "NO\n");
//     Eigen::Matrix3d A;
//     A << 3, 2, 1, 2, 3, 1, 1, 1, 3;

//     std::cout << "\tIS  symmetric: "
//               << (systemVarMatrix_.isApprox(systemVarMatrix_.transpose())
//                       ? "YES\n"
//                       : "NO\n");

//     for (int i = 0; i < myMesh_.numAllNodes_; i++)
//     {
//         // use kdtree find indexes of neighbor nodes

//         for (int j = 0; j < myMesh_.numAllNodes_; j++)
//         {
//             if (systemVarMatrix_.coeffRef(i, j) !=
//                 systemVarMatrix_.coeffRef(j, i))
//             {
//                 std::cout << myMesh_.getNodes()[i][0] << ", "
//                           << myMesh_.getNodes()[i][1] << " -> "
//                           << systemVarMatrix_.coeffRef(i, j) << ", "
//                           << systemVarMatrix_.coeffRef(j, i) << std::endl;
//             }
//         }
//     }

//     // std::cout << systemVarMatrix_ << std::endl;
// }

// void SimulationDomain::solveDomain()
// {
//     Eigen::VectorXd rhs = Eigen::VectorXd::Zero(myMesh_.numAllNodes_);
//     Eigen::VectorXd x0(myMesh_.numAllNodes_);

//     solution_.resize(myMesh_.numAllNodes_);

//     for (int i = 0; i < myMesh_.numAllNodes_; i++)
//     {
//         if (i >= myMesh_.numInnNodes_ && i < myMesh_.numInnNodes_ + 30)
//         {
//             rhs(i) = 100;
//         }
//         x0(i) = 100.0;
//     }

//     // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>>
//     // solver;

//     // solver.compute(systemVarMatrix_);

//     // if (solver.info() != Eigen::Success)
//     // {
//     //     std::cout << " decomposition failed" << std::endl;
//     //     // return;
//     // }
//     // // solution_ = solver.solveWithGuess(rhs, x0);
//     // solution_ = solver.solve(rhs);

//     // std::cout << "#iterations:     " << solver.iterations() << std::endl;
//     // std::cout << "estimated error: " << solver.error() << std::endl;

//     // solution_ = solver.solve(rhs);

//     Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
//         solver;

//     // Compute the ordering permutation vector from the structural pattern
//     // of A
//     solver.analyzePattern(systemVarMatrix_);
//     // Compute the numerical factorization
//     solver.factorize(systemVarMatrix_);
//     // Use the factors to solve the linear system
//     solution_ = solver.solve(rhs_);
// }

// void SimulationDomain::exportData()
// {
//     char fileoutput[256] = "output.txt";

//     std::ofstream myfout(fileoutput);

//     std::vector<double> VX(2);

//     for (int i = 0; i < myMesh_.numAllNodes_; i++)
//     {
//         VX = myMesh_.getNodes()[i];

//         myfout << VX[0] << " " << VX[1] << " " << solution_(i) << std::endl;
//     }

//     myfout.close();
// }

// // explicit instantiation, put this at end of file
// template class SimulationDomain<Rectangle, MQBasis2D>;