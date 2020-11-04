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
