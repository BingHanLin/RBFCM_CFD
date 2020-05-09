#include "simulationDomain.hpp"
#include "MQBasis.hpp"
#include "enumMap.hpp"
#include "pugixml.hpp"
#include "rectangle.hpp"
#include "vtkFileIO.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

SimulationDomain::SimulationDomain(std::shared_ptr<MeshData> mesh,
                                   std::shared_ptr<MQBasis> RBFBasis)
    : controlData_(controlData::instance()),
      myMesh_(mesh),
      myRBFBasis_(RBFBasis),
      viscous_(0.0),
      density_(0.0),
      tStepSize_(0.0),
      endTime_(0.0),
      theta_(0.0),
      convectionVel_(),
      diffusionCoeff_(0.0),
      neighborNum_(0)
{
    setUpSimulation();
    showSummary();
    setupLinearSystem();
    assembleCoeffMatrix();
    assembleRhs();
    solveDomain();
    writeDataToVTK();
}

void SimulationDomain::setupLinearSystem()
{
    std::cout << "#setupLinearSystem" << std::endl;

    varCoeffMatrix_.resize(myMesh_->numOfNodes(), myMesh_->numOfNodes());
    varRhs_.resize(myMesh_->numOfNodes());
    preVarSol_.resize(myMesh_->numOfNodes());

    laplaceMatrix_.resize(myMesh_->numOfNodes(), myMesh_->numOfNodes());
    dxMatrix_.resize(myMesh_->numOfNodes(), myMesh_->numOfNodes());
    dyMatrix_.resize(myMesh_->numOfNodes(), myMesh_->numOfNodes());
    dzMatrix_.resize(myMesh_->numOfNodes(), myMesh_->numOfNodes());

    for (int nodeID = 0; nodeID < myMesh_->numOfNodes(); ++nodeID)
    {
        auto cloud = myMesh_->neighborNodesCloud(nodeID, neighborNum_);

        Eigen::VectorXd laplaceVector =
            myRBFBasis_->collectOnNodes(cloud.nodes, rbfOperatorType::LAPLACE);
        Eigen::VectorXd dxVector = myRBFBasis_->collectOnNodes(
            cloud.nodes, rbfOperatorType::PARTIAL_D1);
        Eigen::VectorXd dyVector = myRBFBasis_->collectOnNodes(
            cloud.nodes, rbfOperatorType::PARTIAL_D2);
        Eigen::VectorXd dzVector = myRBFBasis_->collectOnNodes(
            cloud.nodes, rbfOperatorType::PARTIAL_D3);

        for (int i = 0; i < neighborNum_; i++)
        {
            laplaceMatrix_.insert(nodeID, cloud.id[i]) = laplaceVector[i];
            dxMatrix_.insert(nodeID, cloud.id[i]) = dxVector[i];
            dyMatrix_.insert(nodeID, cloud.id[i]) = dyVector[i];
            dzMatrix_.insert(nodeID, cloud.id[i]) = dzVector[i];
        }
    }
}

void SimulationDomain::setUpSimulation()
{
    std::cout << "#setUpSimulation" << std::endl;

    const auto solverControls = controlData_->paramsDataAt({"solverConstrol"});
    neighborNum_ = solverControls.at("neighborNumber");
    tStepSize_ = solverControls.at("timeStepSize");
    endTime_ = solverControls.at("endTime");
    solverType_ = solverControls.at("solverType");
    theta_ = solverControls.at("transferEqOptions").at("theta");

    const auto physicalControls =
        controlData_->paramsDataAt({"physicsControl"});
    diffusionCoeff_ =
        physicalControls.at("transferEqOptions").at("diffusionCoeff");
    convectionVel_ =
        physicalControls.at("transferEqOptions").at("convectionVel");
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
}

void SimulationDomain::assembleCoeffMatrix()
{
    std::cout << "#assembleCoeffMatrix" << std::endl;

    varCoeffMatrix_.data().squeeze();
    varCoeffMatrix_.reserve(
        Eigen::VectorXi::Constant(myMesh_->numOfNodes(), neighborNum_));

    for (int nodeID = 0; nodeID < myMesh_->numOfNodes(); ++nodeID)
    {
        auto cloud = myMesh_->neighborNodesCloud(nodeID, neighborNum_);

        if (myMesh_->nodeBC(nodeID) == nullptr)
        {
            Eigen::VectorXd localVector = myRBFBasis_->collectOnNodes(
                cloud.nodes, rbfOperatorType::CONSTANT);

            localVector += (-theta_ * tStepSize_ * diffusionCoeff_ *
                            myRBFBasis_->collectOnNodes(
                                cloud.nodes, rbfOperatorType::LAPLACE));

            for (int i = 0; i < neighborNum_; i++)
            {
                varCoeffMatrix_.insert(nodeID, cloud.id[i]) = localVector[i];
            }
        }
        else
        {
            myMesh_->nodeBC(nodeID)->fillCoeffMatrix(nodeID, cloud, myRBFBasis_,
                                                     varCoeffMatrix_);
        }
    }

    std::cout << "#assembleCoeffMatrix end" << std::endl;
}

void SimulationDomain::assembleRhs()
{
    std::cout << "#assembleRhs" << std::endl;

    Eigen::VectorXd rhsInnerVector =
        ((1 - theta_) * tStepSize_ * diffusionCoeff_ * laplaceMatrix_ -
         tStepSize_ * convectionVel_[0] * dxMatrix_ -
         tStepSize_ * convectionVel_[1] * dyMatrix_ -
         tStepSize_ * convectionVel_[2] * dzMatrix_) *
        preVarSol_;

    for (int nodeID = 0; nodeID < myMesh_->numOfNodes(); ++nodeID)
    {
        auto cloud = myMesh_->neighborNodesCloud(nodeID, neighborNum_);

        if (myMesh_->nodeBC(nodeID) == nullptr)
        {
            varRhs_(nodeID) = rhsInnerVector(nodeID);
        }
        else
        {
            myMesh_->nodeBC(nodeID)->fillRhsVector(nodeID, cloud, myRBFBasis_,
                                                   varRhs_);
        }
    }

    std::cout << "#assembleRhs end" << std::endl;
}

void SimulationDomain::solveDomain()
{
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

    // varSol_ = solver.solve(rhs);

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
        solver;

    // Compute the ordering permutation vector from the structural
    // pattern of A
    solver.analyzePattern(varCoeffMatrix_);
    // Compute the numerical factorization
    solver.factorize(varCoeffMatrix_);
    // Use the factors to solve the linear system
    varSol_ = solver.solve(varRhs_);

    std::cout << "#solveDomain" << std::endl;
}

void SimulationDomain::writeDataToVTK() const
{
    pugi::xml_document doc;
    pugi::xml_node VTKFile = doc.append_child("VTKFile");
    VTKFile.append_attribute("type") = "UnstructuredGrid";
    VTKFile.append_attribute("version") = "0.1";
    VTKFile.append_attribute("byte_order") = "LittleEndian";
    pugi::xml_node UnstructuredGrid = VTKFile.append_child("UnstructuredGrid");

    pugi::xml_node Piece = UnstructuredGrid.append_child("Piece");
    Piece.append_attribute("NumberOfPoints") = myMesh_->numOfNodes();
    Piece.append_attribute("NumberOfCells") = myMesh_->numOfNodes();

    pugi::xml_node Points = Piece.append_child("Points");
    appendArrayToVTKNode(myMesh_->nodes(), "Position", Points);

    pugi::xml_node PointData = Piece.append_child("PointData");
    appendScalarsToVTKNode(varSol_, "Variable", PointData);

    pugi::xml_node Cells = Piece.append_child("Cells");
    addCells(myMesh_->numOfNodes(), Cells);

    std::filesystem::create_directories(controlData_->vtkDir());
    const std::string fileName =
        controlData_->vtkDir().string() + "/" + "test" + ".vtu";
    doc.save_file(fileName.c_str());
}