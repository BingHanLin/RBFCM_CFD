#ifndef VTKFILEIO_HPP
#define VTKFILEIO_HPP

#include "pugixml.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <fstream>
template <typename T, int N>
void appendArrayToVTKNode(const std::vector<Eigen::Matrix<T, N, 1>> &values,
                          const std::string name, pugi::xml_node &node)
{
    pugi::xml_node dataArray = node.append_child("DataArray");
    std::stringstream positions;
    for (size_t i = 0; i < values.size(); ++i)
    {
        for (size_t j = 0; j < N; ++j)
        {
            positions << values[i][j] << " ";
        }
    }
    dataArray.append_child(pugi::node_pcdata)
        .set_value(positions.str().c_str());
    dataArray.append_attribute("Name") = name.c_str();

    if (std::is_same<T, double>::value)
    {
        dataArray.append_attribute("type") = "Float64";
    }
    else if (std::is_same<T, size_t>::value)
    {
        dataArray.append_attribute("type") = "Int64";
    }

    dataArray.append_attribute("NumberOfComponents") =
        std::to_string(N).c_str();
    dataArray.append_attribute("format") = "ascii";
}

template <typename T>
inline void appendScalarsToVTKNode(
    const Eigen::Matrix<T, Eigen::Dynamic, 1> &values, const std::string name,
    pugi::xml_node &node)
{
    pugi::xml_node dataArray = node.append_child("DataArray");
    std::stringstream valuesStream;
    for (size_t i = 0; i < values.size(); ++i)
    {
        valuesStream << values[i] << " ";
    }
    dataArray.append_child(pugi::node_pcdata)
        .set_value(valuesStream.str().c_str());
    dataArray.append_attribute("Name") = name.c_str();

    if (std::is_same<T, double>::value)
    {
        dataArray.append_attribute("type") = "Float64";
    }
    else if (std::is_same<T, size_t>::value)
    {
        dataArray.append_attribute("type") = "Int64";
    }

    dataArray.append_attribute("format") = "ascii";
}

void addCells(const size_t cellNum, pugi::xml_node &Node);

void writeVTKGroupFile(const std::string gourpFileName,
                       const std::string childFileNmae, const double timeStep);
#endif
