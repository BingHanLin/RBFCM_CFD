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

    if (paramsStream.good()) paramsStream >> paramsData_;
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
