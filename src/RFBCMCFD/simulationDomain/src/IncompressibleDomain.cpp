#include "IncompressibleDomain.hpp"
#include "ControlData.hpp"
#include "InitialCondition.hpp"
#include "MQBasis.hpp"
#include "MeshData.hpp"
#include "PUConditionPool.hpp"
#include "constant.hpp"
#include "enumMap.hpp"
#include "freeFunctions.hpp"

#include "Rectangle.hpp"
#include "pugixml.hpp"
#include "vtkFileIO.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

IncompressibleDomain::IncompressibleDomain(ControlData* controlData,
                                           MQBasis* RBFBasis,
                                           MeshData* meshData,
                                           PUConditionPool* conditionPool)
    : controlData_(controlData),
      RBFBasis_(RBFBasis),
      meshData_(meshData),
      conditionPool_(conditionPool)
{
    setupSimulation();
    showSummary();
    clearVTKDirectory();
    setupLinearSystem();
    initializeField();
}

void IncompressibleDomain::setupSimulation()
{
    std::cout << "#setupSimulation" << std::endl;

    const auto solverControls = controlData_->paramsDataAt({"solverControl"});

    neighborRadius_ = solverControls.at("neighborRadius");
    tStepSize_ = solverControls.at("timeStepSize");
    endTime_ = solverControls.at("endTime");
    writeInterval_ = solverControls.at("writeInterval");
    dim_ = solverControls.at("dimension");

    const auto physicalControls =
        controlData_->paramsDataAt({"physicsControl"});
    viscosity_ = physicalControls.at("viscosity");
    density_ = physicalControls.at("density");
}

void IncompressibleDomain::showSummary()
{
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
    std::cout << std::setfill(' ') << std::setw(50) << "Summary of simulation"
              << std::endl;
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::setfill(' ')
              << std::endl;

    std::cout << "Number of nodes: " << std::setw(8) << meshData_->numOfNodes()
              << std::endl;

    std::cout << "Time step size: " << std::setw(8) << tStepSize_ << std::endl;
    std::cout << "End time: " << std::setw(8) << endTime_ << std::endl;
    std::cout << "Write Interval: " << std::setw(8) << writeInterval_
              << std::endl;
    std::cout << "Neighbor radius: " << std::setw(8) << neighborRadius_
              << std::endl;

    std::cout << "Viscosity: " << std::setw(8) << viscosity_ << std::endl;
    std::cout << "Density: " << std::setw(8) << density_ << std::endl;
}

void IncompressibleDomain::setupLinearSystem()
{
    std::cout << "#setupLinearSystem" << std::endl;

    const auto numOfNodes = meshData_->numOfNodes();

    // numOfNodes x numOfNodes
    laplaceMatrix_.data().squeeze();
    laplaceMatrix_.reserve(
        Eigen::VectorXi::Constant(numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        auto cloud = meshData_->cloudByID(nodeID);
        auto nodes = meshData_->nodes();

        Eigen::VectorXd laplaceVector =
            RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::LAPLACE);

        for (size_t i = 0; i < cloud.size_; i++)
        {
            laplaceMatrix_.insert(nodeID, cloud.ids_[i]) = laplaceVector[i];
        }
    }

    // dim_* numOfNodes x numOfNodes
    firstOrderDerMatrix_.data().squeeze();
    firstOrderDerMatrix_.reserve(
        Eigen::VectorXi::Constant(dim_ * numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    for (size_t d = 0; d < dim_; ++d)
    {
        const auto start = d * numOfNodes;

        for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
        {
            const auto cloud = meshData_->cloudByID(nodeID);

            Eigen::VectorXd firstDerVector = RBFBasis_->collectOnNodes(
                start + nodeID, firstOrderOperatorTypes[d]);

            for (size_t i = 0; i < cloud.size_; i++)
            {
                firstOrderDerMatrix_.insert(start + nodeID, cloud.ids_[i]) =
                    firstDerVector[i];
            }
        }
    }
}

void IncompressibleDomain::initializeField()
{
    const auto numOfNodes = meshData_->numOfNodes();

    velSol_.resize(dim_ * numOfNodes);
    pSol_.resize(numOfNodes);

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->PIC())
        {
            conditionPool_->PIC()->fillVector(nodeID, pSol_);
        }

        if (conditionPool_->UIC())
        {
            conditionPool_->UIC()->fillVector(nodeID, velSol_);
        }
    }
}

void IncompressibleDomain::assembleCoeffMatrix()
{
    std::cout << "#assembleCoeffMatrix" << std::endl;

    const auto numOfNodes = meshData_->numOfNodes();

    velCoeffMatrix_.data().squeeze();
    velCoeffMatrix_.reserve(
        Eigen::VectorXi::Constant(numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->UBCByNodeID(nodeID) != nullptr)
        {
            conditionPool_->UBCByNodeID(nodeID)->fillCoeffMatrix(
                nodeID, RBFBasis_, velCoeffMatrix_);
        }
        else
        {
            velCoeffMatrix_.row(nodeID) = laplaceMatrix_.row(nodeID);

            velCoeffMatrix_.coeffRef(nodeID, nodeID) -=
                2.0 / viscosity_ / tStepSize_;
        }
    }

    phiCoeffMatrix_.data().squeeze();
    phiCoeffMatrix_.reserve(
        Eigen::VectorXi::Constant(numOfNodes, ESTIMATE_NEIGHBOR_NUM));
    bool refPhiGiven = false;

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->UBCByNodeID(nodeID) != nullptr)
        {
            const auto cloud = meshData_->cloudByID(nodeID);

            Eigen::VectorXd localVector =
                RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::NEUMANN);
            for (size_t i = 0; i < localVector.size(); i++)
            {
                phiCoeffMatrix_.insert(nodeID, cloud.ids_[i]) = localVector[i];
            }
        }
        else
        {
            if (refPhiGiven == false)
            {
                phiCoeffMatrix_.insert(nodeID, nodeID) = 1.0;
                refPhiGiven = true;
            }
        }
    }
}

void IncompressibleDomain::assembleRhs()
{
    // std::cout << "#assembleRhs" << std::endl;

    // Eigen::VectorXd rhsInnerVector =
    //     Eigen::VectorXd::Zero(meshData_->numOfNodes());

    // if (systemSateType_ == systemSateType::TRANSIENT)
    // {
    //     rhsInnerVector = preVarSol_;

    //     rhsInnerVector +=
    //         ((1 - theta_) * tStepSize_ * diffusionCoeff_ * laplaceMatrix_) *
    //         preVarSol_;

    //     rhsInnerVector +=
    //         ((1 - theta_) * tStepSize_ * convectionVel_[0] * dxMatrix_ +
    //          (1 - theta_) * tStepSize_ * convectionVel_[1] * dyMatrix_ +
    //          (1 - theta_) * tStepSize_ * convectionVel_[2] * dzMatrix_) *
    //         preVarSol_;
    // }

    // for (size_t nodeID = 0; nodeID < meshData_->numOfNodes(); ++nodeID)
    // {
    //     auto cloud = meshData_->cloudByID(nodeID);

    //     if (conditionPool_->BCByNodeID(nodeID) == nullptr)
    //     {
    //         varRhs_(nodeID) = rhsInnerVector(nodeID);
    //     }
    //     else
    //     {
    //         conditionPool_->BCByNodeID(nodeID)->fillRhsVector(nodeID,
    //         RBFBasis_,
    //                                                           varRhs_);
    //     }
    // }
}

void IncompressibleDomain::solveMatrix()
{
    // std::cout << "#solveMatrix" << std::endl;

    // // =========================================
    // // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double,
    // Eigen::RowMajor>> solver;
    // // solver.compute(varCoeffMatrix_);

    // // if (solver.info() != Eigen::Success)
    // // {
    // //     std::cout << " decomposition failed" << std::endl;
    // //     return;
    // // }
    // // // solution_ = solver.solveWithGuess(rhs, x0);
    // // varSol_ = solver.solve(varRhs_);

    // // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // // std::cout << "estimated error: " << solver.error() << std::endl;

    // // varSol_ = solver.solve(varRhs_);

    // // =========================================
    // Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>,
    // Eigen::COLAMDOrdering<int>>
    //     solver;

    // // Compute the ordering permutation vector from the structural
    // // pattern of A
    // solver.analyzePattern(varCoeffMatrix_);
    // // Compute the numerical factorization
    // solver.factorize(varCoeffMatrix_);
    // // Use the factors to solve the linear system
    // varSol_ = solver.solve(varRhs_);
    // // =========================================

    // // fill A and b
    // // Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::RowMajor>,
    // //                          Eigen::Lower | Eigen::Upper>
    // //     cg;
    // // cg.compute(varCoeffMatrix_);
    // // varSol_ = cg.solve(varRhs_);
    // // std::cout << "#iterations:     " << cg.iterations() << std::endl;
    // // std::cout << "estimated error: " << cg.error() << std::endl;
    // // // update b, and solve again
    // // varSol_ = cg.solve(varRhs_);
    // // =========================================

    // preVarSol_ = varSol_;
}

Eigen::VectorXd IncompressibleDomain::crankNicolsonU(
    const Eigen::VectorXd& prePSol, const Eigen::VectorXd& preVelSol)
{
    Eigen::VectorXd temp1(preVelSol);
    Eigen::VectorXd temp2(preVelSol);
    Eigen::VectorXd temp3(preVelSol);

    double currError;
    size_t iterNum = 0;

    do
    {
        temp2 = innerCrankNicolsonU(prePSol, preVelSol, temp1);

        temp3 = temp1 - temp2;

        currError = temp3.lpNorm<1>() / temp2.lpNorm<1>();

        temp1 = temp2;

        std::cout << "CrankNicolson error at " << iterNum << " = " << std::fixed
                  << std::setprecision(10) << currError << std::endl;

        iterNum++;

        if (iterNum > crankNicolsonMaxIter_) break;

    } while (currError > crankNicolsonEpsilon_);

    return temp2;
}

Eigen::VectorXd IncompressibleDomain::innerCrankNicolsonU(
    const Eigen::VectorXd& prePSol, const Eigen::VectorXd& preVelSol,
    const Eigen::VectorXd& velTemp)
{
    Eigen::VectorXd velHalf = velTemp * 0.5 + preVelSol * 0.5;

    Eigen::VectorXd RHS = Eigen::VectorXd::Zero(preVelSol.size());

    const auto numOfNodes = meshData_->numOfNodes();
    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->UBCByNodeID(nodeID) != nullptr)
        {
            conditionPool_->UBCByNodeID(nodeID)->fillRhsVector(nodeID,
                                                               RBFBasis_, RHS);
        }
        else
        {
            for (int d = 0; d < dim_; ++d)
            {
                const auto start = d * numOfNodes;
                const auto end = (d + 1) * numOfNodes;

                for (int dd = 0; dd < dim_; ++dd)
                {
                    double velGrad =
                        firstOrderDerMatrix_.row(dd * numOfNodes + nodeID) *
                        velHalf.segment(start, end);

                    RHS(start + nodeID) +=
                        velGrad * velHalf(dd * numOfNodes + nodeID);
                }

                RHS(start + nodeID) =
                    2.0 * RHS(start + nodeID) / viscosity_ -
                    (2.0 / tStepSize_ / viscosity_) * velHalf(start + nodeID);

                double pGrad = firstOrderDerMatrix_.row(start + nodeID) *
                               prePSol.segment(start, end);
                RHS(start + nodeID) += pGrad * (2.0 / viscosity_);

                double velLaplacian = laplaceMatrix_.row(start + nodeID) *
                                      preVelSol.segment(start, end);
                RHS(start + nodeID) -= velLaplacian;
            }
        }
    }

    auto velOut = Eigen::VectorXd::Zero(preVelSol.size());

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                    Eigen::COLAMDOrdering<int>>
        solver;

    // Compute the ordering permutation vector from the structural
    // pattern of A
    solver.analyzePattern(velCoeffMatrix_);
    // Compute the numerical factorization
    solver.factorize(velCoeffMatrix_);
    // Use the factors to solve the linear system
    for (int d = 0; d < dim_; ++d)
    {
        const auto start = d * numOfNodes;
        const auto end = (d + 1) * numOfNodes;
        velOut.segment(start, end) = solver.solve(RHS.segment(start, end));
    }

    return velOut;
}

void IncompressibleDomain::solveDomain()
{
    assembleCoeffMatrix();

    writeDataToVTK();

    auto preVelSol = velSol_;
    auto prePSol = pSol_;

    while (currentTime_ < endTime_)
    {
        crankNicolsonU(prePSol, preVelSol);

        assembleRhs();
        solveMatrix();
        currentTime_ += tStepSize_;
        if (remainder(currentTime_, writeInterval_) <= 0) writeDataToVTK();
    }
}

void IncompressibleDomain::clearVTKDirectory() const
{
    if (std::filesystem::exists(controlData_->vtkDir()))
    {
        for (const auto& entry :
             std::filesystem::directory_iterator(controlData_->vtkDir()))
            std::filesystem::remove_all(entry.path());
    }
    else
    {
        std::filesystem::create_directories(controlData_->vtkDir());
    }
}

void IncompressibleDomain::writeDataToVTK() const
{
    // pugi::xml_document doc;
    // pugi::xml_node VTKFile = doc.append_child("VTKFile");
    // VTKFile.append_attribute("type") = "UnstructuredGrid";
    // VTKFile.append_attribute("version") = "0.1";
    // VTKFile.append_attribute("byte_order") = "LittleEndian";
    // pugi::xml_node UnstructuredGrid =
    // VTKFile.append_child("UnstructuredGrid");

    // pugi::xml_node Piece = UnstructuredGrid.append_child("Piece");
    // Piece.append_attribute("NumberOfPoints") = meshData_->numOfNodes();
    // Piece.append_attribute("NumberOfCells") = meshData_->numOfNodes();

    // pugi::xml_node Points = Piece.append_child("Points");
    // appendArrayToVTKNode(meshData_->nodes(), "Position", Points);

    // pugi::xml_node PointData = Piece.append_child("PointData");
    // appendScalarsToVTKNode(varSol_, "Variable", PointData);

    // pugi::xml_node Cells = Piece.append_child("Cells");
    // addCells(meshData_->numOfNodes(), Cells);

    // // std::filesystem::create_directories(vtkDir());

    // const std::string childFileNmae = controlData_->vtkDir().string() +
    // "/" +
    //                                   std::to_string(currentTime_) +
    //                                   ".vtu";
    // doc.save_file(childFileNmae.c_str());

    // const std::string relChildFileNmae = std::to_string(currentTime_) +
    // ".vtu";

    // const std::string gourpFileName =
    //     controlData_->vtkDir().string() + "/result.pvd";
    // writeVTKGroupFile(gourpFileName, relChildFileNmae, currentTime_);
}
