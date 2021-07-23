#include "controlData.hpp"
#include <fstream>
#include <iostream>

ControlData::ControlData()
{
    setWorkingDir(std::filesystem::current_path());
    readParamsData();
}

ControlData::ControlData(const std::filesystem::path workingDir)
{
    setWorkingDir(workingDir);
    readParamsData();
}

void ControlData::setWorkingDir(const std::filesystem::path workingDir)
{
    workingDir_ = workingDir;
    std::cout << "set working directory: " << workingDir_.string() << '\n';
}

std::filesystem::path ControlData::workingDir() const
{
    return workingDir_;
}

std::filesystem::path ControlData::vtkDir() const
{
    return workingDir().concat("/vtk");
}

void ControlData::readParamsData()
{
    std::ifstream paramsStream(workingDir().concat("/params.json").string());
    if (paramsStream.good())
    {
        paramsStream >> paramsData_;

        const auto solverControls = paramsData_.at("solverControl");
        estimateNeighborNum_ = solverControls.at("estimateNeighborNumber");
        neighborRadius_ = solverControls.at("neighborRadius");
        tStepSize_ = solverControls.at("timeStepSize");
        endTime_ = solverControls.at("endTime");
        writeInterval_ = solverControls.at("writeInterval");
        systemSateType_ = solverControls.at("systemSateType");
        dim_ = solverControls.at("dimension");
        solverType_ = solverControls.at("solverType");
        theta_ = solverControls.at("theta");

        const auto physicalControls = paramsData_.at("physicsControl");
        diffusionCoeff_ = physicalControls.at("diffusionCoeff");
        convectionVel_ = physicalControls.at("convectionVel");
    }
}

const nlohmann::json ControlData::paramsDataAt(
    const std::vector<std::string> searchEntries) const
{
    nlohmann::json queryJSON = paramsData_;

    for (auto oneEntry : searchEntries)
    {
        queryJSON = queryJSON.at(oneEntry);
    }
    return queryJSON;
}

const std::vector<std::string> ControlData::groupNames() const
{
    std::vector<std::string> groupNameList;

    for (const auto& [BCType, BCData] :
         paramsData_.at("boundaryConditions").items())
    {
        for (const auto& oneBCData : BCData)
        {
            groupNameList.push_back(oneBCData.at("groupName"));
        }
    }

    return groupNameList;
}
