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
#include "freeFunctions.hpp"
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
    crankNicolsonEpsilon_ = solverControls.at("crankNicolsonEpsilon");
    crankNicolsonMaxIter_ = solverControls.at("crankNicolsonMaxIter");

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

    std::cout << "CrankNicolson Epsilon: " << std::setw(8)
              << crankNicolsonEpsilon_ << std::endl;

    std::cout << "CrankNicolson Max Iteration: " << std::setw(8)
              << crankNicolsonMaxIter_ << std::endl;

    std::cout << "Viscosity: " << std::setw(8) << viscosity_ << std::endl;
    std::cout << "Density: " << std::setw(8) << density_ << std::endl;
}

void IncompressibleDomain::setupLinearSystem()
{
    std::cout << "#setupLinearSystem" << std::endl;

    const auto numOfNodes = meshData_->numOfNodes();

    // numOfNodes x numOfNodes
    // laplaceMatrix_.data().squeeze();
    // laplaceMatrix_.reserve(
    //     Eigen::VectorXi::Constant(numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    laplaceMatrix_.resize(numOfNodes, numOfNodes);
    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        auto cloud = meshData_->cloudByID(nodeID);

        const Eigen::VectorXd laplaceVector =
            RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::LAPLACE);

        for (size_t i = 0; i < cloud.size_; i++)
        {
            laplaceMatrix_.insert(nodeID, cloud.ids_[i]) = laplaceVector[i];
        }
    }

    // dim_* numOfNodes x numOfNodes firstOrderDerMatrix_.data().squeeze();
    // firstOrderDerMatrix_.reserve(
    //     Eigen::VectorXi::Constant(dim_ * numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    firstOrderDerMatrix_.resize(numOfNodes * dim_, numOfNodes);

    for (size_t d = 0; d < dim_; ++d)
    {
        const auto start = d * numOfNodes;

        for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
        {
            const auto cloud = meshData_->cloudByID(nodeID);

            const Eigen::VectorXd firstDerVector =
                RBFBasis_->collectOnNodes(nodeID, firstOrderOperatorTypes_[d]);

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

    // velCoeffMatrix_.data().squeeze();
    // velCoeffMatrix_.reserve(
    //     Eigen::VectorXi::Constant(numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    velCoeffMatrix_.resize(numOfNodes, numOfNodes);

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

    // phiCoeffMatrix_.data().squeeze();
    // phiCoeffMatrix_.reserve(
    //     Eigen::VectorXi::Constant(numOfNodes, ESTIMATE_NEIGHBOR_NUM));

    phiCoeffMatrix_.resize(numOfNodes, numOfNodes);

    bool refPhiGiven = false;

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->UBCByNodeID(nodeID) != nullptr)
        {
            if (refPhiGiven == false)
            {
                phiCoeffMatrix_.insert(nodeID, nodeID) = 1.0;
                refPhiGiven = true;
            }
            else
            {
                const auto cloud = meshData_->cloudByID(nodeID);

                Eigen::VectorXd localVector =
                    RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::NEUMANN);
                for (size_t i = 0; i < localVector.size(); i++)
                {
                    phiCoeffMatrix_.insert(nodeID, cloud.ids_[i]) =
                        localVector[i];
                }
            }
        }
        else
        {
            phiCoeffMatrix_.row(nodeID) = laplaceMatrix_.row(nodeID);
        }
    }
}

void IncompressibleDomain::solveDomain()
{
    std::cout << "IncompressibleDomain::solveDomain" << std::endl;

    assembleCoeffMatrix();

    writeDataToVTK();

    Eigen::VectorXd velSolNext = velSol_;
    Eigen::VectorXd pSolNext = pSol_;

    while (currentTime_ < endTime_)
    {
        currentTime_ += tStepSize_;

        Eigen::VectorXd velS = crankNicolsonU(pSol_, velSol_);
        solvePU(pSol_, velS, pSolNext, velSolNext);

        velSol_ = velSolNext;
        pSol_ = pSolNext;

        if (writeNow(currentTime_, writeInterval_, tStepSize_))
            writeDataToVTK();
    }
}

Eigen::VectorXd IncompressibleDomain::crankNicolsonU(
    const Eigen::VectorXd& pSol, const Eigen::VectorXd& velSol) const
{
    Eigen::VectorXd temp1(velSol);
    Eigen::VectorXd temp2(velSol);
    Eigen::VectorXd temp3(velSol);

    double currError;
    size_t iterNum = 0;

    do
    {
        temp2 = innerCrankNicolsonU(pSol, velSol, temp1);

        temp3 = temp1 - temp2;

        currError = temp3.lpNorm<1>() / temp2.lpNorm<1>();

        temp1 = temp2;

        iterNum++;

        if (iterNum > crankNicolsonMaxIter_) break;

    } while (currError > crankNicolsonEpsilon_);

    std::cout << "CrankNicolson iterNum: " << iterNum
              << ", error: " << currError << std::endl;

    return temp2;
}

Eigen::VectorXd IncompressibleDomain::innerCrankNicolsonU(
    const Eigen::VectorXd& pSol, const Eigen::VectorXd& velSol,
    const Eigen::VectorXd& velTemp) const
{
    Eigen::VectorXd velHalf = velTemp * 0.5 + velSol * 0.5;

    Eigen::VectorXd RHS = Eigen::VectorXd::Zero(velSol.size());

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

                for (int dd = 0; dd < dim_; ++dd)
                {
                    const double velGrad =
                        firstOrderDerMatrix_.row(dd * numOfNodes + nodeID) *
                        velHalf.segment(start, numOfNodes);

                    RHS(start + nodeID) +=
                        velGrad * velHalf(dd * numOfNodes + nodeID);
                }

                RHS(start + nodeID) =
                    2.0 * RHS(start + nodeID) / viscosity_ -
                    (2.0 / tStepSize_ / viscosity_) * velHalf(start + nodeID);

                const double pGrad =
                    firstOrderDerMatrix_.row(start + nodeID) * pSol;
                RHS(start + nodeID) += pGrad * (2.0 / viscosity_ / density_);

                const double velLaplacian = laplaceMatrix_.row(nodeID) *
                                            velSol.segment(start, numOfNodes);
                RHS(start + nodeID) -= velLaplacian;
            }
        }
    }

    Eigen::VectorXd velOut = Eigen::VectorXd::Zero(velSol.size());

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
        velOut.segment(start, numOfNodes) =
            solver.solve(RHS.segment(start, numOfNodes));
    }

    return velOut;
}

void IncompressibleDomain::solvePU(const Eigen::VectorXd& pSol,
                                   const Eigen::VectorXd& velS,
                                   Eigen::VectorXd& pSolNext,
                                   Eigen::VectorXd& velSolNext) const
{
    const auto numOfNodes = meshData_->numOfNodes();

    pSolNext = pSol;
    velSolNext = velS;

    Eigen::VectorXd RHS = Eigen::VectorXd::Zero(numOfNodes);

    bool refPhiGiven = false;

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->UBCByNodeID(nodeID) != nullptr)
        {
            if (refPhiGiven == false)
            {
                RHS(nodeID) = 0.0;  // reference Phi value
                refPhiGiven = true;
            }
            else
            {
                RHS(nodeID) = 0.0;
            }
        }
        else
        {
            for (int d = 0; d < dim_; ++d)
            {
                const auto start = d * numOfNodes;

                const double velSGrad =
                    firstOrderDerMatrix_.row(start + nodeID) *
                    velS.segment(start, numOfNodes);

                RHS(nodeID) += velSGrad / tStepSize_;
            }
        }
    }

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::RowMajor>,
                    Eigen::COLAMDOrdering<int>>
        solver;

    // Compute the ordering permutation vector from the structural
    // pattern of A
    solver.analyzePattern(phiCoeffMatrix_);
    // Compute the numerical factorization
    solver.factorize(phiCoeffMatrix_);
    // Use the factors to solve the linear system
    const Eigen::VectorXd phi = solver.solve(RHS);

    pSolNext += phi;
    pSolNext -= (viscosity_ * tStepSize_ / 2.0) * laplaceMatrix_ * phi;

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->UBCByNodeID(nodeID) == nullptr)
        {
            for (int d = 0; d < dim_; ++d)
            {
                const auto start = d * numOfNodes;

                const double phiGrad =
                    firstOrderDerMatrix_.row(start + nodeID) * phi;

                velSolNext(start + nodeID) -= phiGrad * tStepSize_ / density_;
            }
        }
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
    const auto numOfNodes = meshData_->numOfNodes();

    pugi::xml_document doc;
    pugi::xml_node VTKFile = doc.append_child("VTKFile");
    VTKFile.append_attribute("type") = "UnstructuredGrid";
    VTKFile.append_attribute("version") = "0.1";
    VTKFile.append_attribute("byte_order") = "LittleEndian";
    pugi::xml_node UnstructuredGrid = VTKFile.append_child("UnstructuredGrid");

    pugi::xml_node Piece = UnstructuredGrid.append_child("Piece");
    Piece.append_attribute("NumberOfPoints") = numOfNodes;
    Piece.append_attribute("NumberOfCells") = numOfNodes;

    pugi::xml_node Points = Piece.append_child("Points");
    appendArrayToVTKNode(meshData_->nodes(), "Position", Points);

    pugi::xml_node PointData = Piece.append_child("PointData");

    const Eigen::VectorXd uSol = velSol_.segment(0, numOfNodes);
    const Eigen::VectorXd vSol = velSol_.segment(numOfNodes, numOfNodes);
    // const Eigen::VectorXd wSol =
    //     velSol_.segment(numOfNodes * 2, numOfNodes * 3);

    appendScalarsToVTKNode(uSol, "u", PointData);
    appendScalarsToVTKNode(vSol, "v", PointData);
    // appendScalarsToVTKNode(wSol, "w", PointData);
    appendScalarsToVTKNode(pSol_, "p", PointData);

    pugi::xml_node Cells = Piece.append_child("Cells");
    addCells(numOfNodes, Cells);

    // std::filesystem::create_directories(vtkDir());

    const std::string childFileNmae = controlData_->vtkDir().string() + "/" +
                                      std::to_string(currentTime_) + ".vtu";
    doc.save_file(childFileNmae.c_str());

    const std::string relChildFileNmae = std::to_string(currentTime_) + ".vtu";

    const std::string gourpFileName =
        controlData_->vtkDir().string() + "/result.pvd";
    writeVTKGroupFile(gourpFileName, relChildFileNmae, currentTime_);
}
