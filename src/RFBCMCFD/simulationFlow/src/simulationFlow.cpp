#include "simulationFlow.hpp"
#include "MQBasis.hpp"
#include "domainData.hpp"
#include "enumMap.hpp"
#include "freeFunctions.hpp"
#include "initialCondition.hpp"
#include "pugixml.hpp"
#include "rectangle.hpp"
#include "vtkFileIO.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

SimulationFlow::SimulationFlow(std::shared_ptr<controlData> inControlData,
                               std::shared_ptr<DomainData> domainData,
                               std::shared_ptr<MQBasis> RBFBasis)
    : controlData_(inControlData),
      myDomainData_(domainData),
      myRBFBasis_(RBFBasis)

{
    setupSimulation();
    showSummary();
    clearVTKDirectory();
    setupLinearSystem();
    initializeField();
}

void SimulationFlow::setupSimulation()
{
    std::cout << "#setupSimulation" << std::endl;
}

void SimulationFlow::showSummary()
{
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::endl;
    std::cout << std::setfill(' ') << std::setw(50) << "Summary of simulation"
              << std::endl;
    std::cout << std::setfill('=') << std::setw(80) << "=" << std::setfill(' ')
              << std::endl;

    std::cout << "Number of nodes: " << std::setw(8)
              << myDomainData_->meshData()->numOfNodes() << std::endl;

    std::cout << "Time step size: " << std::setw(8) << controlData_->tStepSize_
              << std::endl;
    std::cout << "End time: " << std::setw(8) << controlData_->endTime_
              << std::endl;
    std::cout << "Write Interval: " << std::setw(8)
              << controlData_->writeInterval_ << std::endl;
    std::cout << "Neighbor radius: " << std::setw(8)
              << controlData_->neighborRadius_ << std::endl;
    std::cout << "Estimate neighbor number: " << std::setw(8)
              << controlData_->estimateNeighborNum_ << std::endl;

    std::cout << "Diffusity: " << std::setw(8) << controlData_->diffusionCoeff_
              << std::endl;
    std::cout << "Convectivity: " << std::setw(8)
              << controlData_->convectionVel_[0] << ", "
              << controlData_->convectionVel_[1] << ", "
              << controlData_->convectionVel_[2] << std::endl;
}

void SimulationFlow::setupLinearSystem()
{
    std::cout << "#setupLinearSystem" << std::endl;

    const auto numOfNodes = myDomainData_->meshData()->numOfNodes();

    varCoeffMatrix_.resize(numOfNodes, numOfNodes);
    varRhs_.resize(numOfNodes);
    preVarSol_.resize(numOfNodes);

    laplaceMatrix_.resize(numOfNodes, numOfNodes);
    dxMatrix_.resize(numOfNodes, numOfNodes);
    dyMatrix_.resize(numOfNodes, numOfNodes);
    dzMatrix_.resize(numOfNodes, numOfNodes);

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        auto cloud = myDomainData_->meshData()->cloudByID(nodeID);
        auto nodes = myDomainData_->meshData()->nodes();

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

void SimulationFlow::initializeField()
{
    const auto numOfNodes = myDomainData_->meshData()->numOfNodes();

    varSol_ = Eigen::VectorXd::Constant(numOfNodes, 0.0);

    for (size_t nodeID = 0; nodeID < numOfNodes; ++nodeID)
    {
        if (myDomainData_->ICByID(nodeID))
        {
            myDomainData_->ICByID(nodeID)->fillVector(nodeID, varSol_);
        }
    }

    preVarSol_ = varSol_;
}

void SimulationFlow::assembleCoeffMatrix()
{
    std::cout << "#assembleCoeffMatrix" << std::endl;

    varCoeffMatrix_.data().squeeze();
    varCoeffMatrix_.reserve(
        Eigen::VectorXi::Constant(myDomainData_->meshData()->numOfNodes(),
                                  controlData_->estimateNeighborNum_));

    for (size_t nodeID = 0; nodeID < myDomainData_->meshData()->numOfNodes();
         ++nodeID)
    {
        auto cloud = myDomainData_->meshData()->cloudByID(nodeID);
        auto nodes = myDomainData_->meshData()->nodes();

        if (myDomainData_->BCByID(nodeID) == nullptr)
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
            myDomainData_->BCByID(nodeID)->fillCoeffMatrix(nodeID, myRBFBasis_,
                                                           varCoeffMatrix_);
        }
    }
}

void SimulationFlow::assembleRhs()
{
    std::cout << "#assembleRhs" << std::endl;

    Eigen::VectorXd rhsInnerVector =
        Eigen::VectorXd::Zero(myDomainData_->meshData()->numOfNodes());

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

    for (size_t nodeID = 0; nodeID < myDomainData_->meshData()->numOfNodes();
         ++nodeID)
    {
        auto cloud = myDomainData_->meshData()->cloudByID(nodeID);

        if (myDomainData_->BCByID(nodeID) == nullptr)
        {
            varRhs_(nodeID) = rhsInnerVector(nodeID);
        }
        else
        {
            myDomainData_->BCByID(nodeID)->fillRhsVector(nodeID, myRBFBasis_,
                                                         varRhs_);
        }
    }
}

void SimulationFlow::solveMatrix()
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

void SimulationFlow::solveDomain()
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

void SimulationFlow::clearVTKDirectory() const
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

void SimulationFlow::writeDataToVTK() const
{
    pugi::xml_document doc;
    pugi::xml_node VTKFile = doc.append_child("VTKFile");
    VTKFile.append_attribute("type") = "UnstructuredGrid";
    VTKFile.append_attribute("version") = "0.1";
    VTKFile.append_attribute("byte_order") = "LittleEndian";
    pugi::xml_node UnstructuredGrid = VTKFile.append_child("UnstructuredGrid");

    pugi::xml_node Piece = UnstructuredGrid.append_child("Piece");
    Piece.append_attribute("NumberOfPoints") =
        myDomainData_->meshData()->numOfNodes();
    Piece.append_attribute("NumberOfCells") =
        myDomainData_->meshData()->numOfNodes();

    pugi::xml_node Points = Piece.append_child("Points");
    appendArrayToVTKNode(myDomainData_->meshData()->nodes(), "Position",
                         Points);

    pugi::xml_node PointData = Piece.append_child("PointData");
    appendScalarsToVTKNode(varSol_, "Variable", PointData);

    pugi::xml_node Cells = Piece.append_child("Cells");
    addCells(myDomainData_->meshData()->numOfNodes(), Cells);

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
