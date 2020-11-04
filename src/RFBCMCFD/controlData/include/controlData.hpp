#ifndef CONTROLDATA_HPP
#define CONTROLDATA_HPP

#include <filesystem>

#include "json.h"

class controlData
{
   public:
    controlData();
    controlData(const std::filesystem::path workingDir);

    std::filesystem::path workingDir() const;
    std::filesystem::path vtkDir() const;
    nlohmann::json paramsDataAt(const std::vector<std::string> searchEntries);

   private:
    void setWorkingDir(const std::filesystem::path workingDir);
    void readParamsData();

    std::filesystem::path workingDir_;
    nlohmann::json paramsData_;
};

#endif
