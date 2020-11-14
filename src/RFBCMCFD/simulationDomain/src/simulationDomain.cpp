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

SimulationDomain::SimulationDomain(std::shared_ptr<controlData> inControlData,
                                   std::shared_ptr<MeshData> mesh,
                                   std::shared_ptr<MQBasis> RBFBasis)
    : controlData_(inControlData), myMesh_(mesh), myRBFBasis_(RBFBasis)

{
    setupSimulation();
    showSummary();
    clearVTKDirectory();
    setupLinearSystem();
    setupInitialCondition();
}

void SimulationDomain::setupSimulation()
{
    std::cout << "#setupSimulation" << std::endl;
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
    std::cout << "Time step size: " << std::setw(8) << controlData_->tStepSize_
              << std::endl;
    std::cout << "End time: " << std::setw(8) << controlData_->endTime_
              << std::endl;
    std::cout << "Write Interval: " << std::setw(8)
              << controlData_->writeInterval_ << std::endl;
    std::cout << "Neighbor number: " << std::setw(8)
              << controlData_->neighborNum_ << std::endl;

    std::cout << std::endl;

    std::cout << "Diffusity: " << std::setw(8) << controlData_->diffusionCoeff_
              << std::endl;
    std::cout << "Convectivity: " << std::setw(8)
              << controlData_->convectionVel_[0] << ", "
              << controlData_->convectionVel_[1] << ", "
              << controlData_->convectionVel_[2] << std::endl;
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

    for (size_t nodeID = 0; nodeID < myMesh_->numOfNodes(); ++nodeID)
    {
        auto cloud = myMesh_->nodesCloudByID(nodeID);
        auto nodes = myMesh_->nodes();

        Eigen::VectorXd laplaceVector =
            myRBFBasis_->collectOnNodes(cloud, nodes, rbfOperatorType::LAPLACE);
        Eigen::VectorXd dxVector = myRBFBasis_->collectOnNodes(
            cloud, nodes, rbfOperatorType::PARTIAL_D1);
        Eigen::VectorXd dyVector = myRBFBasis_->collectOnNodes(
            cloud, nodes, rbfOperatorType::PARTIAL_D2);
        Eigen::VectorXd dzVector = myRBFBasis_->collectOnNodes(
            cloud, nodes, rbfOperatorType::PARTIAL_D3);

        for (size_t i = 0; i < cloud.size_; i++)
        {
            laplaceMatrix_.insert(nodeID, cloud.ids_[i]) = laplaceVector[i];
            dxMatrix_.insert(nodeID, cloud.ids_[i]) = dxVector[i];
            dyMatrix_.insert(nodeID, cloud.ids_[i]) = dyVector[i];
            dzMatrix_.insert(nodeID, cloud.ids_[i]) = dzVector[i];
        }
    }
}

void SimulationDomain::setupInitialCondition()
{
    const auto initCondition =
        controlData_->paramsDataAt({"physicsControl", "initialConditions"});

    for (auto& oneBCData : initCondition.at("constantValue"))
    {
        const std::string groupName = oneBCData.at("groupName");
        const double val = oneBCData.at("value");
        if (groupName.empty())
        {
            preVarSol_ = Eigen::VectorXd::Constant(myMesh_->numOfNodes(), val);
        }
        else
        {
            for (size_t& nodeID : myMesh_->groupToNodesMap().at(groupName))
            {
                varSol_[nodeID] = val;
            }
        }
    }

    preVarSol_ = varSol_;
}

void SimulationDomain::assembleCoeffMatrix()
{
    std::cout << "#assembleCoeffMatrix" << std::endl;

    varCoeffMatrix_.data().squeeze();
    varCoeffMatrix_.reserve(Eigen::VectorXi::Constant(
        myMesh_->numOfNodes(), controlData_->neighborNum_));

    for (size_t nodeID = 0; nodeID < myMesh_->numOfNodes(); ++nodeID)
    {
        auto cloud = myMesh_->nodesCloudByID(nodeID);
        auto nodes = myMesh_->nodes();

        if (myMesh_->nodeBC(nodeID) == nullptr)
        {
            Eigen::VectorXd localVector;
            if (controlData_->systemSateType_ == systemSateType::STEADY)
            {
                localVector = controlData_->diffusionCoeff_ *
                              myRBFBasis_->collectOnNodes(
                                  cloud, nodes, rbfOperatorType::LAPLACE);

                localVector += controlData_->convectionVel_[0] *
                               myRBFBasis_->collectOnNodes(
                                   cloud, nodes, rbfOperatorType::PARTIAL_D1);

                localVector += controlData_->convectionVel_[1] *
                               myRBFBasis_->collectOnNodes(
                                   cloud, nodes, rbfOperatorType::PARTIAL_D2);

                localVector += controlData_->convectionVel_[2] *
                               myRBFBasis_->collectOnNodes(
                                   cloud, nodes, rbfOperatorType::PARTIAL_D3);
            }
            else
            {
                localVector = myRBFBasis_->collectOnNodes(
                    cloud, nodes, rbfOperatorType::CONSTANT);

                localVector +=
                    (-controlData_->theta_ * controlData_->tStepSize_ *
                     controlData_->diffusionCoeff_ *
                     myRBFBasis_->collectOnNodes(cloud, nodes,
                                                 rbfOperatorType::LAPLACE));

                localVector +=
                    (-controlData_->theta_ * controlData_->tStepSize_ *
                     controlData_->convectionVel_[0] *
                     myRBFBasis_->collectOnNodes(cloud, nodes,
                                                 rbfOperatorType::PARTIAL_D1));

                localVector +=
                    (-controlData_->theta_ * controlData_->tStepSize_ *
                     controlData_->convectionVel_[1] *
                     myRBFBasis_->collectOnNodes(cloud, nodes,
                                                 rbfOperatorType::PARTIAL_D2));

                localVector +=
                    (-controlData_->theta_ * controlData_->tStepSize_ *
                     controlData_->convectionVel_[2] *
                     myRBFBasis_->collectOnNodes(cloud, nodes,
                                                 rbfOperatorType::PARTIAL_D3));
            }

            for (size_t i = 0; i < cloud.size_; i++)
                varCoeffMatrix_.insert(nodeID, cloud.ids_[i]) = localVector[i];
        }
        else
        {
            myMesh_->nodeBC(nodeID)->fillCoeffMatrix(nodeID, myRBFBasis_,
                                                     varCoeffMatrix_);
        }
    }
}

void SimulationDomain::assembleRhs()
{
    std::cout << "#assembleRhs" << std::endl;

    Eigen::VectorXd rhsInnerVector =
        Eigen::VectorXd::Zero(myMesh_->numOfNodes());

    if (controlData_->systemSateType_ == systemSateType::TRANSIENT)
    {
        rhsInnerVector = preVarSol_;

        rhsInnerVector +=
            ((1 - controlData_->theta_) * controlData_->tStepSize_ *
             controlData_->diffusionCoeff_ * laplaceMatrix_) *
            preVarSol_;

        rhsInnerVector +=
            ((1 - controlData_->theta_) * controlData_->tStepSize_ *
                 controlData_->convectionVel_[0] * dxMatrix_ +
             (1 - controlData_->theta_) * controlData_->tStepSize_ *
                 controlData_->convectionVel_[1] * dyMatrix_ +
             (1 - controlData_->theta_) * controlData_->tStepSize_ *
                 controlData_->convectionVel_[2] * dzMatrix_) *
            preVarSol_;
    }

    for (size_t nodeID = 0; nodeID < myMesh_->numOfNodes(); ++nodeID)
    {
        auto cloud = myMesh_->nodesCloudByID(nodeID);

        if (myMesh_->nodeBC(nodeID) == nullptr)
        {
            varRhs_(nodeID) = rhsInnerVector(nodeID);
        }
        else
        {
            myMesh_->nodeBC(nodeID)->fillRhsVector(nodeID, myRBFBasis_,
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

    if (controlData_->systemSateType_ == systemSateType::STEADY)
    {
        assembleRhs();
        solveMatrix();
        writeDataToVTK();
    }
    else
    {
        writeDataToVTK();

        while (controlData_->currentTime_ < controlData_->endTime_)
        {
            assembleRhs();
            solveMatrix();
            controlData_->currentTime_ += controlData_->tStepSize_;
            if (remainder(controlData_->currentTime_,
                          controlData_->writeInterval_) <= 0)
                writeDataToVTK();
        }
    }
}

void SimulationDomain::clearVTKDirectory() const
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

    // std::filesystem::create_directories(controlData_->vtkDir());

    const std::string childFileNmae =
        controlData_->vtkDir().string() + "/" +
        std::to_string(controlData_->currentTime_) + ".vtu";
    doc.save_file(childFileNmae.c_str());

    const std::string relChildFileNmae =
        std::to_string(controlData_->currentTime_) + ".vtu";

    const std::string gourpFileName =
        controlData_->vtkDir().string() + "/result.pvd";
    writeVTKGroupFile(gourpFileName, relChildFileNmae,
                      controlData_->currentTime_);
}
