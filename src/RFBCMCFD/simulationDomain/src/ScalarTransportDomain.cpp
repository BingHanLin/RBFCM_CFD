#include "ScalarTransportDomain.hpp"
#include "ControlData.hpp"
#include "InitialCondition.hpp"
#include "MQBasis.hpp"
#include "MeshData.hpp"
#include "ScalarConditionPool.hpp"
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

ScalarTransportDomain::ScalarTransportDomain(ControlData* controlData,
                                             MQBasis* RBFBasis,
                                             MeshData* meshData,
                                             ScalarConditionPool* conditionPool)
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

void ScalarTransportDomain::setupSimulation()
{
    std::cout << "#setupSimulation" << std::endl;

    const auto solverControls = controlData_->paramsDataAt({"solverControl"});

    neighborRadius_ = solverControls.at("neighborRadius");
    tStepSize_ = solverControls.at("timeStepSize");
    endTime_ = solverControls.at("endTime");
    writeInterval_ = solverControls.at("writeInterval");
    systemSateType_ = solverControls.at("systemSateType");
    dim_ = solverControls.at("dimension");
    theta_ = solverControls.at("theta");

    const auto physicalControls =
        controlData_->paramsDataAt({"physicsControl"});
    diffusionCoeff_ = physicalControls.at("diffusionCoeff");
    convectionVel_ = physicalControls.at("convectionVel");
}

void ScalarTransportDomain::showSummary()
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

    std::cout << "Diffusity: " << std::setw(8) << diffusionCoeff_ << std::endl;
    std::cout << "Convectivity: " << std::setw(8) << convectionVel_[0] << ", "
              << convectionVel_[1] << ", " << convectionVel_[2] << std::endl;
}

void ScalarTransportDomain::setupLinearSystem()
{
    std::cout << "#setupLinearSystem" << std::endl;

    const auto numOfNodes = meshData_->numOfNodes();

    varCoeffMatrix_.resize(numOfNodes, numOfNodes);
    varRhs_.resize(numOfNodes);
    preVarSol_.resize(numOfNodes);

    laplaceMatrix_.resize(numOfNodes, numOfNodes);
    dxMatrix_.resize(numOfNodes, numOfNodes);
    dyMatrix_.resize(numOfNodes, numOfNodes);
    dzMatrix_.resize(numOfNodes, numOfNodes);

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        auto cloud = meshData_->cloudByID(nodeID);
        auto nodes = meshData_->nodes();

        Eigen::VectorXd laplaceVector =
            RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::LAPLACE);
        Eigen::VectorXd dxVector =
            RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::PARTIAL_D1);
        Eigen::VectorXd dyVector =
            RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::PARTIAL_D2);
        Eigen::VectorXd dzVector =
            RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::PARTIAL_D3);

        for (size_t i = 0; i < cloud.size_; i++)
        {
            laplaceMatrix_.insert(nodeID, cloud.ids_[i]) = laplaceVector[i];
            dxMatrix_.insert(nodeID, cloud.ids_[i]) = dxVector[i];
            dyMatrix_.insert(nodeID, cloud.ids_[i]) = dyVector[i];
            dzMatrix_.insert(nodeID, cloud.ids_[i]) = dzVector[i];
        }
    }
}

void ScalarTransportDomain::initializeField()
{
    const auto numOfNodes = meshData_->numOfNodes();

    varSol_ = Eigen::VectorXd::Constant(numOfNodes, 0.0);

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (conditionPool_->ICByNodeID(nodeID))
        {
            conditionPool_->ICByNodeID(nodeID)->fillVector(nodeID, varSol_);
        }
    }

    preVarSol_ = varSol_;
}

void ScalarTransportDomain::assembleCoeffMatrix()
{
    std::cout << "#assembleCoeffMatrix" << std::endl;

    varCoeffMatrix_.data().squeeze();
    varCoeffMatrix_.reserve(Eigen::VectorXi::Constant(meshData_->numOfNodes(),
                                                      ESTIMATENEIGHBORNUM));

    for (size_t nodeID = 0; nodeID < meshData_->numOfNodes(); ++nodeID)
    {
        auto cloud = meshData_->cloudByID(nodeID);
        auto nodes = meshData_->nodes();

        if (conditionPool_->BCByNodeID(nodeID) == nullptr)
        {
            Eigen::VectorXd localVector;
            if (systemSateType_ == systemSateType::STEADY)
            {
                localVector =
                    diffusionCoeff_ *
                    RBFBasis_->collectOnNodes(nodeID, rbfOperatorType::LAPLACE);

                localVector += convectionVel_[0] *
                               RBFBasis_->collectOnNodes(
                                   nodeID, rbfOperatorType::PARTIAL_D1);

                localVector += convectionVel_[1] *
                               RBFBasis_->collectOnNodes(
                                   nodeID, rbfOperatorType::PARTIAL_D2);

                localVector += convectionVel_[2] *
                               RBFBasis_->collectOnNodes(
                                   nodeID, rbfOperatorType::PARTIAL_D3);
            }
            else
            {
                localVector = RBFBasis_->collectOnNodes(
                    nodeID, rbfOperatorType::CONSTANT);

                localVector += (-theta_ * tStepSize_ * diffusionCoeff_ *
                                RBFBasis_->collectOnNodes(
                                    nodeID, rbfOperatorType::LAPLACE));

                localVector += (-theta_ * tStepSize_ * convectionVel_[0] *
                                RBFBasis_->collectOnNodes(
                                    nodeID, rbfOperatorType::PARTIAL_D1));

                localVector += (-theta_ * tStepSize_ * convectionVel_[1] *
                                RBFBasis_->collectOnNodes(
                                    nodeID, rbfOperatorType::PARTIAL_D2));

                localVector += (-theta_ * tStepSize_ * convectionVel_[2] *
                                RBFBasis_->collectOnNodes(
                                    nodeID, rbfOperatorType::PARTIAL_D3));
            }

            for (size_t i = 0; i < cloud.size_; i++)
                varCoeffMatrix_.insert(nodeID, cloud.ids_[i]) = localVector[i];
        }
        else
        {
            conditionPool_->BCByNodeID(nodeID)->fillCoeffMatrix(
                nodeID, RBFBasis_, varCoeffMatrix_);
        }
    }
}

void ScalarTransportDomain::assembleRhs()
{
    std::cout << "#assembleRhs" << std::endl;

    Eigen::VectorXd rhsInnerVector =
        Eigen::VectorXd::Zero(meshData_->numOfNodes());

    if (systemSateType_ == systemSateType::TRANSIENT)
    {
        rhsInnerVector = preVarSol_;

        rhsInnerVector +=
            ((1 - theta_) * tStepSize_ * diffusionCoeff_ * laplaceMatrix_) *
            preVarSol_;

        rhsInnerVector +=
            ((1 - theta_) * tStepSize_ * convectionVel_[0] * dxMatrix_ +
             (1 - theta_) * tStepSize_ * convectionVel_[1] * dyMatrix_ +
             (1 - theta_) * tStepSize_ * convectionVel_[2] * dzMatrix_) *
            preVarSol_;
    }

    for (size_t nodeID = 0; nodeID < meshData_->numOfNodes(); ++nodeID)
    {
        auto cloud = meshData_->cloudByID(nodeID);

        if (conditionPool_->BCByNodeID(nodeID) == nullptr)
        {
            varRhs_(nodeID) = rhsInnerVector(nodeID);
        }
        else
        {
            conditionPool_->BCByNodeID(nodeID)->fillRhsVector(nodeID, RBFBasis_,
                                                              varRhs_);
        }
    }
}

void ScalarTransportDomain::solveMatrix()
{
    std::cout << "#solveMatrix" << std::endl;

    // =========================================
    // Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double>> solver;
    // solver.compute(varCoeffMatrix_);

    // if (solver.info() != Eigen::Success)
    // {
    //     std::cout << " decomposition failed" << std::endl;
    //     return;
    // }
    // // solution_ = solver.solveWithGuess(rhs, x0);
    // varSol_ = solver.solve(varRhs_);

    // std::cout << "#iterations:     " << solver.iterations() << std::endl;
    // std::cout << "estimated error: " << solver.error() << std::endl;

    // varSol_ = solver.solve(varRhs_);

    // =========================================
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
        solver;

    // Compute the ordering permutation vector from the structural
    // pattern of A
    solver.analyzePattern(varCoeffMatrix_);
    // Compute the numerical factorization
    solver.factorize(varCoeffMatrix_);
    // Use the factors to solve the linear system
    varSol_ = solver.solve(varRhs_);
    // =========================================

    // fill A and b
    // Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
    //                          Eigen::Lower | Eigen::Upper>
    //     cg;
    // cg.compute(varCoeffMatrix_);
    // varSol_ = cg.solve(varRhs_);
    // std::cout << "#iterations:     " << cg.iterations() << std::endl;
    // std::cout << "estimated error: " << cg.error() << std::endl;
    // // update b, and solve again
    // varSol_ = cg.solve(varRhs_);
    // =========================================

    preVarSol_ = varSol_;
}

void ScalarTransportDomain::solveDomain()
{
    assembleCoeffMatrix();

    if (systemSateType_ == systemSateType::STEADY)
    {
        assembleRhs();
        solveMatrix();
        writeDataToVTK();
    }
    else
    {
        writeDataToVTK();

        while (currentTime_ < endTime_)
        {
            assembleRhs();
            solveMatrix();
            currentTime_ += tStepSize_;
            if (remainder(currentTime_, writeInterval_) <= 0) writeDataToVTK();
        }
    }
}

void ScalarTransportDomain::clearVTKDirectory() const
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

void ScalarTransportDomain::writeDataToVTK() const
{
    pugi::xml_document doc;
    pugi::xml_node VTKFile = doc.append_child("VTKFile");
    VTKFile.append_attribute("type") = "UnstructuredGrid";
    VTKFile.append_attribute("version") = "0.1";
    VTKFile.append_attribute("byte_order") = "LittleEndian";
    pugi::xml_node UnstructuredGrid = VTKFile.append_child("UnstructuredGrid");

    pugi::xml_node Piece = UnstructuredGrid.append_child("Piece");
    Piece.append_attribute("NumberOfPoints") = meshData_->numOfNodes();
    Piece.append_attribute("NumberOfCells") = meshData_->numOfNodes();

    pugi::xml_node Points = Piece.append_child("Points");
    appendArrayToVTKNode(meshData_->nodes(), "Position", Points);

    pugi::xml_node PointData = Piece.append_child("PointData");
    appendScalarsToVTKNode(varSol_, "Variable", PointData);

    pugi::xml_node Cells = Piece.append_child("Cells");
    addCells(meshData_->numOfNodes(), Cells);

    // std::filesystem::create_directories(vtkDir());

    const std::string childFileNmae = controlData_->vtkDir().string() + "/" +
                                      std::to_string(currentTime_) + ".vtu";
    doc.save_file(childFileNmae.c_str());

    const std::string relChildFileNmae = std::to_string(currentTime_) + ".vtu";

    const std::string gourpFileName =
        controlData_->vtkDir().string() + "/result.pvd";
    writeVTKGroupFile(gourpFileName, relChildFileNmae, currentTime_);
}
