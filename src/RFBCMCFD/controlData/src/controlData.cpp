#include "controlData.hpp"
#include <fstream>
#include <iostream>

controlData::controlData()
{
    setWorkingDir(std::filesystem::current_path());
    readParamsData();
}

controlData::controlData(const std::filesystem::path workingDir)
{
    setWorkingDir(workingDir);
    readParamsData();
}

void controlData::setWorkingDir(const std::filesystem::path workingDir)
{
    workingDir_ = workingDir;
    std::cout << "set working directory: " << workingDir_.string() << '\n';
}

std::filesystem::path controlData::workingDir() const
{
    return workingDir_;
}

std::filesystem::path controlData::vtkDir() const
{
    return workingDir().concat("/vtk");
}

void controlData::readParamsData()
{
    std::ifstream paramsStream(workingDir().concat("/params.json").string());
    if (paramsStream.good())
    {
        paramsStream >> paramsData_;

        const auto solverControls = paramsData_.at("solverControl");
        neighborNum_ = solverControls.at("neighborNumber");
        tStepSize_ = solverControls.at("timeStepSize");
        endTime_ = solverControls.at("endTime");
        writeInterval_ = solverControls.at("writeInterval");
        systemSateType_ = solverControls.at("systemSateType");
        dim_ = solverControls.at("dimension");
        solverType_ = solverControls.at("solverType");
        theta_ = solverControls.at("transferEqOptions").at("theta");

        const auto physicalControls = paramsData_.at("physicsControl");
        diffusionCoeff_ =
            physicalControls.at("transferEqOptions").at("diffusionCoeff");
        convectionVel_ =
            physicalControls.at("transferEqOptions").at("convectionVel");
    }
}

nlohmann::json controlData::paramsDataAt(
    const std::vector<std::string> searchEntries)
{
    nlohmann::json queryJSON = paramsData_;

    for (auto oneEntry : searchEntries)
    {
        queryJSON = queryJSON.at(oneEntry);
    }
    return queryJSON;
}
