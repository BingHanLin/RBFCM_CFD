#ifndef CONTROLDATA_HPP
#define CONTROLDATA_HPP

#include "json.h"
#include <filesystem>

class ControlData
{
   public:
    ControlData();
    ControlData(const std::filesystem::path workingDir);

    std::filesystem::path workingDir() const;
    std::filesystem::path vtkDir() const;

    const nlohmann::json paramsDataAt(
        const std::vector<std::string> searchEntries) const;

   private:
    void setWorkingDir(const std::filesystem::path workingDir);
    void readParamsData();

    std::filesystem::path workingDir_;
    nlohmann::json paramsData_;
};

#endif
