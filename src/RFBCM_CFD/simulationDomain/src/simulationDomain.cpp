#include "simulationDomain.hpp"
#include "MQBasis.hpp"
#include "enumMap.hpp"
#include "freeFunctions.hpp"
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
      currentTime_(0.0),
      writeInterval_(0.0),
      theta_(0.0),
      convectionVel_(),
      diffusionCoeff_(0.0),
      neighborNum_(0)
{
    setupSimulation();
    showSummary();
    clearVTKDirectory();
    setupLinearSystem();
    setupInitialCondition();
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

        auto indexs = sortIndexes(cloud.id);

        for (int i = 0; i < neighborNum_; i++)
        {
            laplaceMatrix_.insert(nodeID, cloud.id[indexs[i]]) =
                laplaceVector[indexs[i]];
            dxMatrix_.insert(nodeID, cloud.id[indexs[i]]) = dxVector[indexs[i]];
            dyMatrix_.insert(nodeID, cloud.id[indexs[i]]) = dyVector[indexs[i]];
            dzMatrix_.insert(nodeID, cloud.id[indexs[i]]) = dzVector[indexs[i]];
        }
    }
}

void SimulationDomain::setupSimulation()
{
    std::cout << "#setupSimulation" << std::endl;

    const auto solverControls = controlData_->paramsDataAt({"solverConstrol"});
    neighborNum_ = solverControls.at("neighborNumber");
    tStepSize_ = solverControls.at("timeStepSize");
    endTime_ = solverControls.at("endTime");
    writeInterval_ = solverControls.at("writeInterval");

    solverType_ = solverControls.at("solverType");
    theta_ = solverControls.at("transferEqOptions").at("theta");

    const auto physicalControls =
        controlData_->paramsDataAt({"physicsControl"});
    diffusionCoeff_ =
        physicalControls.at("transferEqOptions").at("diffusionCoeff");
    convectionVel_ =
        physicalControls.at("transferEqOptions").at("convectionVel");
}

void SimulationDomain::setupInitialCondition()
{
    const auto initCondition =
        controlData_->paramsDataAt({"physicsControl", "initialConditions"});

    if (initCondition.at("type") == initTypeEnum::UNIFORM)
    {
        const double val = initCondition.at("uniform").at("value");
        varSol_ = Eigen::VectorXd::Constant(myMesh_->numOfNodes(), val);
        preVarSol_ = Eigen::VectorXd::Constant(myMesh_->numOfNodes(), val);
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
    std::cout << "Write Interval: " << std::setw(8) << writeInterval_
              << std::endl;
    std::cout << "Neighbor number: " << std::setw(8) << neighborNum_
              << std::endl;

    std::cout << std::endl;

    std::cout << "Diffusity: " << std::setw(8) << diffusionCoeff_ << std::endl;
    std::cout << "Convectivity: " << std::setw(8) << convectionVel_[0] << ", "
              << convectionVel_[1] << ", " << convectionVel_[2] << std::endl;
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
            Eigen::VectorXd localVector;
            //  =
            //     Eigen::VectorXd::Zero(myMesh_->numOfNodes());
            if (tStepSize_ == 0.0)
            {
                localVector = diffusionCoeff_ *
                              myRBFBasis_->collectOnNodes(
                                  cloud.nodes, rbfOperatorType::LAPLACE);

                // localVector += convectionVel_[0] *
                //                myRBFBasis_->collectOnNodes(
                //                    cloud.nodes, rbfOperatorType::PARTIAL_D1);

                // localVector += convectionVel_[1] *
                //                myRBFBasis_->collectOnNodes(
                //                    cloud.nodes, rbfOperatorType::PARTIAL_D2);

                // localVector += convectionVel_[2] *
                //                myRBFBasis_->collectOnNodes(
                //                    cloud.nodes, rbfOperatorType::PARTIAL_D3);
            }
            else
            {
                localVector = myRBFBasis_->collectOnNodes(
                    cloud.nodes, rbfOperatorType::CONSTANT);

                localVector += (-theta_ * tStepSize_ * diffusionCoeff_ *
                                myRBFBasis_->collectOnNodes(
                                    cloud.nodes, rbfOperatorType::LAPLACE));

                // localVector += (-theta_ * tStepSize_ * convectionVel_[0] *
                //                 myRBFBasis_->collectOnNodes(
                //                     cloud.nodes,
                //                     rbfOperatorType::PARTIAL_D1));

                // localVector += (-theta_ * tStepSize_ * convectionVel_[1] *
                //                 myRBFBasis_->collectOnNodes(
                //                     cloud.nodes,
                //                     rbfOperatorType::PARTIAL_D2));

                // localVector += (-theta_ * tStepSize_ * convectionVel_[2] *
                //                 myRBFBasis_->collectOnNodes(
                //                     cloud.nodes,
                //                     rbfOperatorType::PARTIAL_D3));
            }

            auto indexs = sortIndexes(cloud.id);

            for (int i = 0; i < neighborNum_; i++)
                varCoeffMatrix_.insert(nodeID, cloud.id[indexs[i]]) =
                    localVector[indexs[i]];
        }
        else
        {
            myMesh_->nodeBC(nodeID)->fillCoeffMatrix(nodeID, cloud, myRBFBasis_,
                                                     varCoeffMatrix_);
        }
    }
}

void SimulationDomain::assembleRhs()
{
    std::cout << "#assembleRhs" << std::endl;

    Eigen::VectorXd rhsInnerVector =
        Eigen::VectorXd::Zero(myMesh_->numOfNodes());

    if (tStepSize_ != 0.0)
    {
        rhsInnerVector = preVarSol_;

        rhsInnerVector +=
            ((1 - theta_) * tStepSize_ * diffusionCoeff_ * laplaceMatrix_) *
            preVarSol_;

        // rhsInnerVector +=
        //     ((1 - theta_) * tStepSize_ * convectionVel_[0] * dxMatrix_ +
        //      (1 - theta_) * tStepSize_ * convectionVel_[1] * dyMatrix_ +
        //      (1 - theta_) * tStepSize_ * convectionVel_[2] * dzMatrix_) *
        //     preVarSol_;
    }

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
}

void SimulationDomain::solveMatrix()
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

void SimulationDomain::solveDomain()
{
    assembleCoeffMatrix();
    writeDataToVTK();

    if (tStepSize_ == 0)
    {
        assembleRhs();
        solveMatrix();
        writeDataToVTK();
    }
    else
    {
        while (currentTime_ < endTime_)
        {
            assembleRhs();
            solveMatrix();
            currentTime_ += tStepSize_;
            if (remainder(currentTime_, writeInterval_) <= 0) writeDataToVTK();
        }
    }
}

void SimulationDomain::clearVTKDirectory() const
{
    for (const auto& entry :
         std::filesystem::directory_iterator(controlData_->vtkDir()))
        std::filesystem::remove_all(entry.path());
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

    const std::string childFileNmae = controlData_->vtkDir().string() + "/" +
                                      std::to_string(currentTime_) + ".vtu";
    doc.save_file(childFileNmae.c_str());

    const std::string gourpFileName =
        controlData_->vtkDir().string() + "/result.pvd";
    writeVTKGroupFile(gourpFileName, childFileNmae, currentTime_);
}
