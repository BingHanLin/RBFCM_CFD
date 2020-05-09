#include "controlData.hpp"
#include <fstream>
#include <iostream>

controlData* controlData::instance_ = nullptr;

controlData* controlData::instance()
{
    if (!instance_)
    {
        instance_ = new controlData;
        instance_->setWorkingDir();
        instance_->readParamsData();
    }
    return instance_;
}

void controlData::setWorkingDir()
{
    workingDir_ = std::filesystem::current_path();
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
