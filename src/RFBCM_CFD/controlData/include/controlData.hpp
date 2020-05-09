#ifndef CONTROLDATA_HPP
#define CONTROLDATA_HPP

#include <filesystem>

#include "json.h"

class controlData
{
   public:
    static controlData* instance();
    std::filesystem::path workingDir() const;
    std::filesystem::path vtkDir() const;
    nlohmann::json paramsDataAt(const std::vector<std::string> searchEntries);

   private:
    controlData(){};                    // private so that it can  not be called
    controlData(controlData const&){};  // copy constructor is private
    // controlData& operator=(
    //     controlData const&){};  // assignment operator is private

    void setWorkingDir();
    void readParamsData();

    std::filesystem::path workingDir_;
    nlohmann::json paramsData_;

    static controlData* instance_;
};

#endif
