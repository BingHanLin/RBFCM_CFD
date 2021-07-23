#ifndef CONTROLDATA_HPP
#define CONTROLDATA_HPP

#include "enumMap.hpp"
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
    const std::vector<std::string> groupNames() const;

    // physicsControl
    double viscous_;
    double density_;
    double diffusionCoeff_;
    std::array<double, 3> convectionVel_;

    // solverConstrol
    size_t estimateNeighborNum_;
    double neighborRadius_;
    double endTime_;
    double tStepSize_;
    double writeInterval_;
    double theta_;
    solverType solverType_;
    systemSateType systemSateType_;
    size_t dim_;
    double currentTime_;

    // size_t crankNicolsonMaxIter_;
    // double crankNicolsonEpsilon_;

   private:
    void setWorkingDir(const std::filesystem::path workingDir);
    void readParamsData();

    std::filesystem::path workingDir_;
    nlohmann::json paramsData_;
};

#endif
